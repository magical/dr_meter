#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "libavformat/avformat.h"
#include "libavutil/error.h"

#define VERSION "0.1"

#define MAX_CHANNELS 32
#define MAX_FRAGMENTS 32768 // more than 24h

// The length of the window over which the RMS and peak are calculated.
// Specified in milliseconds. Don't change this!
#define FRAGMENT_LENGTH 3000

#define FACTOR8 ((sample)1.0 / (sample)(1 << 7))
#define FACTOR16 ((sample)1.0 / (sample)(1 << 15))
#define FACTOR32 ((sample)1.0 / (sample)(1UL << 31))

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))

typedef double sample;

const char throbbler[4] = "/-\\|";

static int debug = 0;

struct track_info {
	int err;
	sample dr;
	sample peak;
	sample rms;
	int duration;
	char *artist;
	char *album;
        char *tracknumber;
	char *title;
};

/******************************************************************************/

struct stream_context {
	AVFormatContext *format_ctx;
	int stream_index; // the stream we are decoding
	AVPacket real_pkt;
	AVPacket pkt;
	AVFrame *frame;
	AVDictionary *opts;
	enum {
		STATE_UNINITIALIZED,
		STATE_INITIALIZED,
		STATE_OPEN,
		STATE_VALID_PACKET,
		STATE_NEED_FLUSH,
		STATE_CLOSED,
	} state;
};

void sc_init(struct stream_context *self) {
	self->format_ctx = NULL;
	self->stream_index = 0;
	self->state = STATE_INITIALIZED;
	self->opts = NULL;
}

int sc_open(struct stream_context *self, const char *filename) {
	int err;

	sc_init(self);
	err = avformat_open_input(&self->format_ctx, filename, NULL, NULL);
	if (err < 0) { return err; }
	err = avformat_find_stream_info(self->format_ctx, NULL);
	if (err < 0) { return err; }

	self->state = STATE_OPEN;

	return 0;
}

/* return the AVCodecContext for the active stream */
AVCodecContext *sc_get_codec(struct stream_context *self) {
	return self->format_ctx->streams[self->stream_index]->codec;
}

void sc_close(struct stream_context *self) {
	if (STATE_OPEN <= self->state && self->state != STATE_CLOSED) {
		av_dict_free(&self->opts);
		avcodec_close(sc_get_codec(self));
		avformat_close_input(&self->format_ctx);
		av_free(self->frame);
		self->state = STATE_CLOSED;
	}
}

bool sc_eof(struct stream_context *self) {
	return self->state == STATE_CLOSED;
}

int sc_start_stream(struct stream_context *self, int stream_index) {
	self->stream_index = stream_index;
	AVCodecContext *ctx = sc_get_codec(self);
	AVCodec *codec = avcodec_find_decoder(ctx->codec_id);
	if (codec == NULL) {
		return AVERROR_DECODER_NOT_FOUND;
	}

	// Allocate a frame.
	self->frame = av_frame_alloc();
	if (self->frame == NULL) {
		return AVERROR(ENOMEM);
	}

	if (codec->id == CODEC_ID_AC3) {
		av_dict_set(&self->opts, "drc_scale", "0", 0);
	}

	/* XXX check codec */
	return avcodec_open2(ctx, codec, &self->opts);
}

// Decode one frame of audio
int sc_get_next_frame(struct stream_context *self) {
	int err;

	if (self->state == STATE_CLOSED) {
		return AVERROR_EOF;
	}

	// Grab a new packet, if necessary
	while (self->state == STATE_OPEN) {
		err = av_read_frame(self->format_ctx, &self->real_pkt);
		if (err == AVERROR_EOF) {
			av_init_packet(&self->pkt);
			self->pkt.data = NULL;
			self->pkt.size = 0;
			self->state = STATE_NEED_FLUSH;
		} else if (err < 0) {
			return err;
		} else if (self->real_pkt.stream_index == self->stream_index) {
			self->pkt = self->real_pkt;
			self->state = STATE_VALID_PACKET;
		} else {
			// we don't care about this frame; try another
			av_packet_unref(&self->real_pkt);
		}
	}

	AVCodecContext *codec_ctx = sc_get_codec(self);

	// Decode the audio.

	// The return value is the number of bytes read from the packet.
	// The codec is not required to read the entire packet, so we may need
	// to keep it around for a while.
	int got_frame;
	//avcodec_get_frame_defaults(self->frame);
	av_frame_unref(self->frame);
	self->frame = av_frame_alloc();

	err = avcodec_decode_audio4(codec_ctx, self->frame, &got_frame, &self->pkt);
	if (err < 0) { return err; }
	int bytes_used = err;

	if (self->state == STATE_VALID_PACKET) {
		if (0 < bytes_used && bytes_used < self->pkt.size) {
			if (debug) {
				printf("partial frame\n");
			}
			self->pkt.data += bytes_used;
			self->pkt.size -= bytes_used;
		} else  {
			self->state = STATE_OPEN;
			av_packet_unref(&self->real_pkt);
		}
	} else if (self->state == STATE_NEED_FLUSH) {
		avcodec_close(codec_ctx);
		self->state = STATE_CLOSED;
	}

	return 0;
}

void **sc_get_buf(struct stream_context *self)
{
	uint8_t **tmp = self->frame->extended_data; // for type checking
	return (void **)tmp;
}

size_t sc_get_samples(struct stream_context *self)
{
	return self->frame->nb_samples;
}

char *sc_get_metadata(struct stream_context *self, const char *key)
{
	AVDictionaryEntry *entry = av_dict_get(self->format_ctx->metadata, key, NULL, 0);
	if (entry != NULL) {
		return entry->value;
	}
	return NULL;
}

/******************************************************************************/

struct dr_meter {
	sample *rms_values[MAX_CHANNELS];
	sample *peak_values[MAX_CHANNELS];

	sample sum[MAX_CHANNELS];
	sample peak[MAX_CHANNELS];

	int channels;
	int sample_rate;
	int sample_fmt;
	int sample_size;

	unsigned long long total_samples;

	size_t fragment; // The index of the current fragment
	size_t fragment_size; // The size of a fragment in samples
	size_t fragment_read; // The number of samples scanned so far
	bool fragment_started;
};

int compare_samples(const void *s1, const void *s2) {
	sample rms1 = *(sample *)s1;
	sample rms2 = *(sample *)s2;
	if (rms1 > rms2) return -1;
	else if (rms1 < rms2) return 1;
	return 0;
}

sample get_sample(void *buf, size_t i, enum AVSampleFormat sample_fmt) {
	switch(sample_fmt) {
	case AV_SAMPLE_FMT_U8:
	case AV_SAMPLE_FMT_U8P:
		return (sample)(((uint8_t *)buf)[i] - 0x80) * FACTOR8;
	case AV_SAMPLE_FMT_S16:
	case AV_SAMPLE_FMT_S16P:
		return (sample)(((int16_t *)buf)[i]) * FACTOR16;
	case AV_SAMPLE_FMT_S32:
	case AV_SAMPLE_FMT_S32P:
		return (sample)(((int32_t *)buf)[i]) * FACTOR32;
	case AV_SAMPLE_FMT_FLT:
	case AV_SAMPLE_FMT_FLTP:
		return (sample)(((float *)buf)[i]);
	case AV_SAMPLE_FMT_DBL:
	case AV_SAMPLE_FMT_DBLP:
		return (sample)(((double *)buf)[i]);
	default:
		return 0.0;
	}
}

sample to_db(const sample linear) {
	return 20.0 * log10(linear);
}

void meter_init(struct dr_meter *self) {
	for (int ch = 0; ch < MAX_CHANNELS; ch++) {
		self->rms_values[ch] = NULL;
		self->peak_values[ch] = NULL;
		self->sum[ch] = 0;
		self->peak[ch] = 0;
	}
	self->channels = 0;
	self->sample_rate = 0;
	self->sample_size = 0;
	self->sample_fmt = 0;
	self->total_samples = 0;
	self->fragment = 0;
	self->fragment_size = 0;
	self->fragment_read = 0;
	self->fragment_started = 0;
}

int meter_start(struct dr_meter *self, int channels, int sample_rate, int sample_fmt) {
	if (channels > MAX_CHANNELS) {
		fprintf(stderr, "FATAL ERROR: Too many channels! Max channels %is.\n", MAX_CHANNELS);
		return 240; // ???
	}

	if (sample_fmt != AV_SAMPLE_FMT_U8 &&
	    sample_fmt != AV_SAMPLE_FMT_U8P &&
	    sample_fmt != AV_SAMPLE_FMT_S16 &&
	    sample_fmt != AV_SAMPLE_FMT_S16P &&
	    sample_fmt != AV_SAMPLE_FMT_S32 &&
	    sample_fmt != AV_SAMPLE_FMT_S32P &&
	    sample_fmt != AV_SAMPLE_FMT_FLT &&
	    sample_fmt != AV_SAMPLE_FMT_FLTP &&
	    sample_fmt != AV_SAMPLE_FMT_DBL &&
	    sample_fmt != AV_SAMPLE_FMT_DBLP) {
		fprintf(stderr, "FATAL ERROR: Unsupported sample format: %s\n", av_get_sample_fmt_name(sample_fmt));
		return 240;
	}
        if (debug) {
            fprintf(stderr, "DEBUG: Sample format: %s\n", av_get_sample_fmt_name(sample_fmt));
        }

	// Allocate RMS and peak storage
	for (int ch = 0; ch < channels; ch++) {
		self->rms_values[ch] =
			malloc(MAX_FRAGMENTS * sizeof(*self->rms_values[ch]));
		self->peak_values[ch] =
			malloc(MAX_FRAGMENTS * sizeof(*self->peak_values[ch]));
		if (self->rms_values[ch] == NULL ||
		    self->peak_values[ch] == NULL) {
			return AVERROR(ENOMEM);
		}
	}
	self->channels = channels;
	self->sample_rate = sample_rate;
	self->sample_fmt = sample_fmt;
	self->sample_size = av_get_bytes_per_sample(sample_fmt);
	if (sample_rate == 44100) {
		self->fragment_size = ((long)sample_rate + 60) * FRAGMENT_LENGTH / 1000;
	} else {
		self->fragment_size = (long)sample_rate * FRAGMENT_LENGTH / 1000;
	}
	assert(self->fragment_size > 0);
	return 0;
}

static int meter_fragment_start(struct dr_meter *self) {
	if (self->fragment >= MAX_FRAGMENTS) {
		fprintf(stderr, "FATAL ERROR: Input too long! Max length %is.\n", MAX_FRAGMENTS*3);
		return 240;
	}

	for (int ch = 0; ch < self->channels; ch++) {
		self->sum[ch] = 0;
		self->peak[ch] = 0;
	}

	self->fragment_read = 0;
	self->fragment_started = true;

	return 0;
}

static void meter_fragment_finish(struct dr_meter *self) {
	for (int ch = 0; ch < self->channels; ch++) {
		self->rms_values[ch][self->fragment] =
			sqrt(2.0 * self->sum[ch] / self->fragment_read);
		self->peak_values[ch][self->fragment] = self->peak[ch];
	}
	self->fragment++;
	self->fragment_started = false;
}

static void dump(void *buf, size_t len, int sample_fmt) {
	size_t i;
	for (i = 0; i < len; i++) {
		fprintf(stdout, "[%zd]: %.17f\n", i, get_sample(buf, i, sample_fmt));
	}
}

static inline void meter_scan_internal(struct dr_meter *self, void **buf, size_t start, size_t end, size_t samples, enum AVSampleFormat sample_fmt) {
	if (buf == NULL) {
		return;
	}
	if (av_sample_fmt_is_planar(sample_fmt)) {
		for (int ch = 0; ch < self->channels; ch++)
		for (size_t i = start; i < end; i++) {
			if (buf[ch] == NULL) {
				printf("missing channel %d data\n", ch);
				exit(1);
			}
			sample value = get_sample(buf[ch], i, sample_fmt);
			//if (isnan(value)) {
			//	value = 0;
			//}
			if (debug && isnan(value)) {
				fprintf(stdout, "%f at [%d][%zd] (%zd)\n", value, ch, i, i + self->fragment_read);
				dump(buf[ch], samples, sample_fmt);
				exit(1);
			}
			self->sum[ch] += value * value;

			value = fabs(value);
			if (self->peak[ch] < value) {
				self->peak[ch] = value;
			}
		}
	} else {
		for (size_t i = start; i < end; i++)
		for (int ch = 0; ch < self->channels; ch++) {
			sample value = get_sample(buf[0], i * self->channels + ch, sample_fmt);
			self->sum[ch] += value * value;

			value = fabs(value);
			if (self->peak[ch] < value) {
				self->peak[ch] = value;
			}
		}
	}
}

/* Feed the meter. Scan a single frame of audio. */
int meter_feed(struct dr_meter *self, void **buf, size_t samples) {
	size_t start = 0, end = samples;
	int err;

	if (debug && 0) {
		fprintf(stdout, "%zd\n", samples);
	}

	while (start < samples) {
		if (!self->fragment_started) {
			err = meter_fragment_start(self);
			if (err) return err;
		}

		size_t fragment_left = self->fragment_size - self->fragment_read;
		end = min(fragment_left, samples);
		#define CASE(fmt) case fmt: meter_scan_internal(self, buf, start, end, samples, fmt); break
		switch (self->sample_fmt) {
		CASE(AV_SAMPLE_FMT_U8);
		CASE(AV_SAMPLE_FMT_U8P);
		CASE(AV_SAMPLE_FMT_S16);
		CASE(AV_SAMPLE_FMT_S16P);
		CASE(AV_SAMPLE_FMT_S32);
		CASE(AV_SAMPLE_FMT_S32P);
		CASE(AV_SAMPLE_FMT_FLT);
		CASE(AV_SAMPLE_FMT_FLTP);
		CASE(AV_SAMPLE_FMT_DBL);
		CASE(AV_SAMPLE_FMT_DBLP);
		default:
			meter_scan_internal(self, buf, start, end, samples, self->sample_fmt);
		}
		#undef CASE
		self->fragment_read += end - start;
		start = end;

		if (self->fragment_size <= self->fragment_read) {
			meter_fragment_finish(self);
		}
	}

	//printf("%lld += %zu\n", self->total_samples, samples);
	self->total_samples += samples;

	return 0;
}

int meter_finish(struct dr_meter *self, struct track_info *info) {
	if (self->fragment_started) {
		// process any leftover audio
		meter_fragment_finish(self);
	}
	sample rms_score[MAX_CHANNELS];
	sample rms[MAX_CHANNELS];
	sample peak_score[MAX_CHANNELS];
	sample dr_channel[MAX_CHANNELS];
	sample dr_sum = 0;
	sample peak_max = 0;
	sample rms_sum = 0;
	size_t values_to_use = 0;
	values_to_use = max(self->fragment / 5, 1);
	if (!values_to_use) {
		info->peak = 0;
		info->rms = 0;
		info->dr = -INFINITY;
		info->duration = 0;
		return 0;
	}
	for (int ch = 0; ch < self->channels; ch++) {
		qsort(self->rms_values[ch], self->fragment, sizeof(*self->rms_values[ch]), compare_samples);
		sample rms_ch_sum = 0;
		for (size_t i = 0; i < values_to_use; i++) {
			sample value = self->rms_values[ch][i];
			rms_ch_sum += value * value;
		}
		rms_score[ch] = sqrt(rms_ch_sum / values_to_use);

		for (size_t i = values_to_use; i < self->fragment; i++) {
			sample value = self->rms_values[ch][i];
			rms_ch_sum += value * value;
		}
		rms[ch] = sqrt(rms_ch_sum / self->fragment);
		rms_sum += rms_ch_sum / self->fragment;

		qsort(self->peak_values[ch], self->fragment, sizeof(*self->peak_values[ch]), compare_samples);
		peak_score[ch] = self->peak_values[ch][min(1, self->fragment)];
		peak_max = max(peak_max, self->peak_values[ch][0]);

		dr_channel[ch] = to_db(peak_score[ch] / rms_score[ch]);
		dr_sum += dr_channel[ch];
		printf("Ch. %2i:  Peak %8.2f (%8.2f)   RMS %8.2f (%8.2f)   DR = %6.2f\n",
		       ch,
		       to_db(self->peak_values[ch][0]),
		       to_db(peak_score[ch]),
		       to_db(rms[ch]),
		       to_db(rms_score[ch]),
		       dr_channel[ch]);
	}
	info->peak = peak_max;
	info->rms = sqrt(rms_sum / (sample)self->channels);
	info->dr = dr_sum / (sample)self->channels;
	info->duration = (self->total_samples + self->sample_rate/2) / self->sample_rate;
	//printf("%lld\n", self->total_samples);
	printf("DR%d\n", (int)round(info->dr));
	return 0;
}

void meter_free(struct dr_meter *self) {
	for (int ch = 0; ch < self->channels; ch++) {
		free(self->rms_values[ch]);
		free(self->peak_values[ch]);
		self->rms_values[ch] = NULL;
		self->peak_values[ch] = NULL;
	}
}

/******************************************************************************/


int print_av_error(const char *function_name, int error) {
	char errorbuf[128];
	char *error_ptr = errorbuf;
	if (av_strerror(error, errorbuf, sizeof(errorbuf)) < 0) {
		error_ptr = strerror(AVUNERROR(error));
	}
	fprintf(stderr, "dr_meter: %s: %s\n", function_name, error_ptr);
	return error;
}

char *strdup(const char *s)
{
	if (s == NULL || *s == '\0')
		return NULL;
	size_t len = strlen(s);
	char *s2 = malloc(len + 1);
	if (s2 == NULL)
		return NULL;
	strcpy(s2, s);
	return s2;
}

int calculate_track_dr(const char *filename, struct track_info *t, int number) {
	struct stream_context sc;
	struct dr_meter meter;
	int err;

	meter_init(&meter);

	err = sc_open(&sc, filename);
	if (err < 0) { return print_av_error("sc_open", err); }

	if (debug) {
		printf("DEBUG: Codec: %s\n", avcodec_get_name(sc_get_codec(&sc)->codec_id));
	}


	int stream_index = err = av_find_best_stream(
		sc.format_ctx, AVMEDIA_TYPE_AUDIO, -1, -1, NULL, 0);
	if (err < 0) { print_av_error("av_find_best_stream", err); goto cleanup; }

	err = sc_start_stream(&sc, stream_index);
	if (err < 0) { print_av_error("sc_start_stream", err); goto cleanup; }

	t->artist = strdup(sc_get_metadata(&sc, "artist"));
	t->album = strdup(sc_get_metadata(&sc, "album"));
        t->tracknumber = strdup(sc_get_metadata(&sc, "track"));
	t->title = strdup(sc_get_metadata(&sc, "title"));

	AVCodecContext *codec_ctx = sc_get_codec(&sc);
	err = meter_start(&meter, codec_ctx->channels, codec_ctx->sample_rate, codec_ctx->sample_fmt);
	if (err) { goto cleanup; }

	size_t fragment = 0;
	int throbbler_stage = 0;
	while (!sc_eof(&sc)) {
		err = sc_get_next_frame(&sc);
		if (err < 0) {
			print_av_error("sc_get_next_frame", err);
			goto cleanup;
		}

		err = meter_feed(&meter, sc_get_buf(&sc), sc_get_samples(&sc));
		if (err) { goto cleanup; }

		if (fragment < meter.fragment) {
			fragment = meter.fragment;
			if ((throbbler_stage % 4) == 0) {
				fprintf(stderr, "\033[1K\033[1G"
						"Analyzing track %i... "
						"%c  %2zu:%02zu ",
						number,
						throbbler[throbbler_stage / 4],
						(fragment * 3) / 60,
						(fragment * 3) % 60);
			}
			throbbler_stage += 1;
			throbbler_stage %= 16;
		}
	}

	fprintf(stderr, "\033[1K\033[1G");
	meter_finish(&meter, t);

cleanup:
	meter_free(&meter);
	sc_close(&sc);

	if (err < 0) {
		return err;
	}

	return 0;
}

void print_bar(int ch)
{
	if (ch == '=') {
		printf("================================================================================\n");
	} else {
		printf("--------------------------------------------------------------------------------\n");
	}
}

bool streq(const char *a, const char *b)
{
	return strcmp(a, b) == 0;
}

// array is a pointer to an array of structures which contain a char* member at the given offset
const char *collapse(void *array, int count, size_t size, size_t offset, const char *zero, const char *many)
{
	const char *ret = zero;
	bool first = true;
	int i;
	for (i = 0; i < count; i++) {
		// s = array[i]->member
		char *s = *(char **)((char *)array + size*i + offset);
		if (s == NULL || s[0] == '\0') {
			continue;
		}

		if (first) {
			ret = s;
			first = false;
		} else if (!streq(ret, s)) {
			ret = many;
			break;
		}
	}
	return ret;
}

void print_dr(struct track_info *info, int count)
{
	int i;
	const char *artist = collapse(info, count, sizeof(struct track_info), offsetof(struct track_info, artist), "Unknown", "Unknown");
	const char *album = collapse(info, count, sizeof(struct track_info), offsetof(struct track_info, album), "Unknown", "Various Artists");
	sample dr_sum = 0;
	printf("dr_meter " VERSION "\n");
	//printf("log date: \n");
	printf("\n");
	print_bar('-');
	printf("Analyzed: %s / %s\n", artist, album);
	print_bar('-');
	printf("\n");
	printf("DR         Peak         RMS     Duration Track\n");
	print_bar('-');
	for (i = 0; i < count; i++) {
		// XXX don't print tracks with errors
		struct track_info *t = &info[i];
		if (t->err) continue;
		printf("DR%-4d %8.2f dB %8.2f dB %6d:%02d %s-%s\n",
		       (int)round(t->dr), to_db(t->peak), to_db(t->rms),
		       t->duration/60, t->duration%60,
		       t->tracknumber, t->title);
		dr_sum += t->dr;
	}
	print_bar('-');
	printf("\n");
	printf("Number of tracks:  %i\n", count);
	printf("Official DR value: DR%i\n", (int)round(dr_sum / count));
	print_bar('=');
	printf("\n");
}

int main(int argc, char** argv) {
	av_register_all();
	av_log_set_level(AV_LOG_ERROR);

	bool err_occurred = false;
	int err;

	if (argc <= 1) {
		fprintf(stderr, "Reading from standard input...\n");
		err = 0;//do_calculate_dr("pipe:");
		if (err) {
			err_occurred = true;
		}
	} else {
		struct track_info *info = calloc(sizeof(struct track_info), argc-1);
		if (info == NULL) {
			exit(1);
		}
		for (int i = 1; i < argc; i++) {
			info->err = calculate_track_dr(argv[i], &info[i-1], i);
			if (info->err) {
				err_occurred = true;
			}
		}
		print_dr(info, argc - 1);
	}

	exit(err_occurred ? EXIT_FAILURE : EXIT_SUCCESS);
}
