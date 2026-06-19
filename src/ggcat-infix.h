#pragma once

#include <zlib.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <vector>

#include "ketopt.h"
#include "kmer.h"
#include "ggcat-colormap.h"

struct GgcatInfixParams {
    int k;
    int n_thread;
};

static inline uint8_t ggcat_comp_base(uint8_t b) {
    return 3 - b;
}

static inline uint64_t ggcat_mask_bits(int len) {
    return len == 32 ? UINT64_MAX : ((1ULL << (2 * len)) - 1ULL);
}

static inline uint64_t ggcat_revcomp_bits(uint64_t x, int len) {
    uint64_t y = 0;
    for (int i = 0; i < len; ++i) {
        y = (y << 2) | (uint64_t)(3 - (x & 3ULL));
        x >>= 2;
    }
    return y;
}

static inline uint64_t ggcat_encode_base(const char *seq, int pos) {
    int c = seq_nt4_table[(uint8_t)seq[pos]];
    if (c >= 4) throw std::runtime_error("ggcat unitig contains a non-ACGT base");
    return (uint64_t)c;
}

struct GgcatCanonicalEdge {
    uint64_t infix;
    uint8_t edge;
};

static inline bool ggcat_canonical_edge(uint64_t kp1mer, int k, GgcatCanonicalEdge *out) {
    uint64_t k_mask = ggcat_mask_bits(k);
    uint64_t left_kmer = kp1mer >> 2;
    uint64_t right_kmer = kp1mer & k_mask;
    if (left_kmer == right_kmer) return false;

    uint64_t rc = ggcat_revcomp_bits(kp1mer, k + 1);
    if (kp1mer == rc) return false;

    uint64_t infix_f = (kp1mer >> 2) & ggcat_mask_bits(k - 1);
    uint64_t infix_r = (rc >> 2) & ggcat_mask_bits(k - 1);
    uint64_t y = infix_f < infix_r ? kp1mer : rc;

    out->infix = (y >> 2) & ggcat_mask_bits(k - 1);
    out->edge = (uint8_t)((((y >> (2 * k)) & 3ULL) << 2) | (y & 3ULL));
    return true;
}

struct GgcatEndpointRecord {
    uint64_t infix;
    uint32_t sid;
    uint8_t side;
    uint8_t base;
};

static inline void ggcat_add_unitig_boundary_endpoints(std::vector<GgcatEndpointRecord>& records,
                                                       const char *seq,
                                                       int len,
                                                       int k,
                                                       uint32_t first_sid,
                                                       uint32_t last_sid) {
    int infix_len = k - 1;
    uint64_t infix_mask = ggcat_mask_bits(infix_len);
    uint64_t first_infix = 0;
    uint8_t first_base = 0;
    for (int i = 0; i < infix_len; ++i) {
        uint64_t c = ggcat_encode_base(seq, i);
        if (i == 0) first_base = (uint8_t)c;
        first_infix = (first_infix << 2) | c;
    }

    uint8_t right_base = (uint8_t)ggcat_encode_base(seq, infix_len);
    uint64_t rc = ggcat_revcomp_bits(first_infix, infix_len);
    if (first_infix < rc) {
        records.push_back(GgcatEndpointRecord{first_infix, first_sid, 1, right_base});
    } else {
        records.push_back(GgcatEndpointRecord{rc, first_sid, 0, ggcat_comp_base(right_base)});
    }

    int left_pos = len - k;
    int last_infix_pos = left_pos + 1;
    uint8_t left_base = left_pos == 0 ? first_base : (uint8_t)ggcat_encode_base(seq, left_pos);
    uint64_t last_infix = first_infix;
    if (last_infix_pos > 0) {
        if (last_infix_pos <= infix_len / 2) {
            for (int i = 0; i < last_infix_pos; ++i) {
                int next_pos = infix_len + i;
                uint64_t c = next_pos == infix_len ? right_base : ggcat_encode_base(seq, next_pos);
                last_infix = ((last_infix << 2) & infix_mask) | c;
            }
        } else {
            last_infix = 0;
            for (int i = 0; i < infix_len; ++i) {
                last_infix = (last_infix << 2) | ggcat_encode_base(seq, last_infix_pos + i);
            }
        }
    }

    rc = ggcat_revcomp_bits(last_infix, infix_len);
    if (last_infix < rc) {
        records.push_back(GgcatEndpointRecord{last_infix, last_sid, 0, left_base});
    } else {
        records.push_back(GgcatEndpointRecord{rc, last_sid, 1, ggcat_comp_base(left_base)});
    }
}

static inline uint64_t ggcat_hist_idx(uint64_t j, uint64_t sigma) {
    return (sigma * (sigma - 1)) / 2 + j - 1;
}

static inline bool ggcat_is_hex(uint8_t c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
}

static inline uint32_t ggcat_hex_value(uint8_t c) {
    if (c <= '9') return (uint32_t)(c - '0');
    if (c <= 'F') return (uint32_t)(c - 'A' + 10);
    return (uint32_t)(c - 'a' + 10);
}

static void ggcat_parse_color_runs_append(const char *s,
                                          size_t n,
                                          std::vector<std::pair<uint32_t, uint32_t> >& runs) {
    size_t pos = 0;
    while (pos + 2 < n) {
        if (s[pos] != 'C' || s[pos + 1] != ':') {
            ++pos;
            continue;
        }
        pos += 2;

        uint32_t sid = 0;
        size_t sid_start = pos;
        while (pos < n && ggcat_is_hex((uint8_t)s[pos])) {
            sid = (sid << 4) | ggcat_hex_value((uint8_t)s[pos]);
            ++pos;
        }
        if (pos == sid_start || pos >= n || s[pos] != ':') continue;
        ++pos;

        uint32_t count = 0;
        size_t count_start = pos;
        while (pos < n && s[pos] >= '0' && s[pos] <= '9') {
            count = 10 * count + (uint32_t)(s[pos] - '0');
            ++pos;
        }
        if (pos == count_start) continue;
        runs.push_back(std::make_pair(sid, count));
    }
}

static inline bool ggcat_endpoint_record_less(const GgcatEndpointRecord& a,
                                              const GgcatEndpointRecord& b) {
    if (a.infix != b.infix) return a.infix < b.infix;
    if (a.side != b.side) return a.side < b.side;
    if (a.base != b.base) return a.base < b.base;
    return a.sid < b.sid;
}

static inline bool ggcat_endpoint_record_equal(const GgcatEndpointRecord& a,
                                               const GgcatEndpointRecord& b) {
    return a.infix == b.infix && a.side == b.side && a.base == b.base && a.sid == b.sid;
}

struct GgcatEndpointScratch {
    explicit GgcatEndpointScratch(size_t n_colors)
        : words((n_colors + 63) / 64),
          seen(words),
          multi(words) {
        for (int i = 0; i < 4; ++i) {
            left[i].resize(words);
            right[i].resize(words);
        }
        for (int i = 0; i < 16; ++i) edge[i].resize(words);
    }

    void clear() {
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(multi.begin(), multi.end(), 0);
        for (int i = 0; i < 4; ++i) {
            std::fill(left[i].begin(), left[i].end(), 0);
            std::fill(right[i].begin(), right[i].end(), 0);
        }
        for (int i = 0; i < 16; ++i) std::fill(edge[i].begin(), edge[i].end(), 0);
    }

    size_t words;
    std::array<std::vector<uint64_t>, 4> left;
    std::array<std::vector<uint64_t>, 4> right;
    std::array<std::vector<uint64_t>, 16> edge;
    std::vector<uint64_t> seen;
    std::vector<uint64_t> multi;
};

static inline void ggcat_or_colors(std::vector<uint64_t>& bits,
                                   const std::vector<uint32_t>& colors) {
    for (size_t i = 0; i < colors.size(); ++i) {
        uint32_t c = colors[i];
        bits[c >> 6] |= 1ULL << (c & 63);
    }
}

static inline uint64_t ggcat_popcount_bits(const std::vector<uint64_t>& bits) {
    uint64_t n = 0;
    for (size_t i = 0; i < bits.size(); ++i) n += (uint64_t)__builtin_popcountll(bits[i]);
    return n;
}

static void ggcat_add_endpoint_group_to_hist(const std::vector<GgcatEndpointRecord>& records,
                                             size_t begin,
                                             size_t end,
                                             int k,
                                             GgcatColorCache& colors,
                                             GgcatEndpointScratch& scratch,
                                             std::vector<uint64_t>& hist) {
    uint64_t infix = records[begin].infix;
    scratch.clear();

    for (size_t i = begin; i < end; ++i) {
        const std::vector<uint32_t>& cs = colors.colors(records[i].sid);
        if (records[i].side == 0) ggcat_or_colors(scratch.left[records[i].base], cs);
        else ggcat_or_colors(scratch.right[records[i].base], cs);
    }

    for (uint8_t a = 0; a < 4; ++a) {
        for (uint8_t b = 0; b < 4; ++b) {
            uint64_t kp1mer = ((uint64_t)a << (2 * k)) | (infix << 2) | (uint64_t)b;
            GgcatCanonicalEdge ce;
            if (!ggcat_canonical_edge(kp1mer, k, &ce)) continue;
            if (ce.infix != infix || ce.edge != (uint8_t)((a << 2) | b)) continue;

            std::vector<uint64_t>& e = scratch.edge[(a << 2) | b];
            for (size_t w = 0; w < scratch.words; ++w) {
                uint64_t x = scratch.left[a][w] & scratch.right[b][w];
                e[w] = x;
                scratch.multi[w] |= scratch.seen[w] & x;
                scratch.seen[w] |= x;
            }
        }
    }

    uint64_t sigma = ggcat_popcount_bits(scratch.seen);
    if (sigma == 0) return;

    for (int e = 0; e < 16; ++e) {
        uint64_t j = 0;
        for (size_t w = 0; w < scratch.words; ++w) {
            j += (uint64_t)__builtin_popcountll(scratch.edge[e][w] & ~scratch.multi[w]);
        }
        if (j > 0) ++hist[ggcat_hist_idx(j, sigma)];
    }
}

static void output_hist_infix_ggcat(int argc, char *argv[]) {
#ifndef PANGROWTH_WITH_GGCAT
    (void)argc;
    (void)argv;
    fprintf(stderr, "pangrowth hist_infix_ggcat requires a build with -DPANGROWTH_WITH_GGCAT=ON\n");
    return;
#else
    GgcatInfixParams param;
    param.k = 17;
    param.n_thread = 1;

    ketopt_t options = KETOPT_INIT;
    int c;
    while ((c = ketopt(&options, argc, argv, 1, "k:t:", 0)) >= 0) {
        if (c == 'k') param.k = atoi(options.arg);
        else if (c == 't') param.n_thread = atoi(options.arg);
    }

    if (argc - options.ind != 1) {
        fprintf(stderr, "Usage: pangrowth hist_infix_ggcat [options] <ggcat_k_graph.fa>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -k INT     k-mer size used to build the ggcat graph [%d]\n", param.k);
        fprintf(stderr, "  -t INT     threads for ggcat colormap queries [%d]\n", param.n_thread);
        return;
    }
    if (param.k < 2 || param.k > 31) {
        fprintf(stderr, "Error: hist_infix_ggcat currently requires 2 <= k <= 31\n");
        exit(EXIT_FAILURE);
    }

    const char *graph_path = argv[options.ind];
    GgcatColorCache color_cache(graph_path, (size_t)param.n_thread);
    uint64_t n_colors = color_cache.num_colors();
    if (n_colors == 0) {
        fprintf(stderr, "Error: ggcat graph has no colors\n");
        exit(EXIT_FAILURE);
    }
    size_t n_color_subsets = color_cache.num_color_subsets();
    if (n_color_subsets == 0) {
        fprintf(stderr, "Error: ggcat graph has no color subsets\n");
        exit(EXIT_FAILURE);
    }

    std::vector<uint64_t> hist((n_colors + 1) * n_colors / 2, 0);
    std::vector<GgcatEndpointRecord> endpoint_records;
    std::vector<uint32_t> subset_ids(n_color_subsets);
    for (size_t i = 0; i < n_color_subsets; ++i) subset_ids[i] = (uint32_t)i;
    color_cache.reserve(n_color_subsets);
    color_cache.query(subset_ids.data(), subset_ids.size());

    gzFile fp = gzopen(graph_path, "r");
    if (fp == 0) {
        fprintf(stderr, "Could not open ggcat graph: %s\n", graph_path);
        exit(EXIT_FAILURE);
    }
    kseq_t *ks = kseq_init(fp);
    int ret;
    uint64_t unitigs = 0, internal_edges = 0;
    uint64_t color_run_count = 0, mixed_internal_edges = 0;
    std::vector<std::pair<uint32_t, uint32_t> > runs;
    while ((ret = kseq_read(ks)) >= 0) {
        int len = (int)ks->seq.l;
        if (len < param.k) continue;
        size_t n_kmers = (size_t)(len - param.k + 1);

        runs.clear();
        if (ks->name.l) ggcat_parse_color_runs_append(ks->name.s, ks->name.l, runs);
        if (ks->comment.l) ggcat_parse_color_runs_append(ks->comment.s, ks->comment.l, runs);
        if (runs.empty()) {
            throw std::runtime_error("ggcat unitig has no color runs; build the graph with ggcat build -c");
        }
        size_t run_kmers = 0;
        for (size_t i = 0; i < runs.size(); ++i) {
            run_kmers += runs[i].second;
        }
        if (run_kmers != n_kmers) {
            throw std::runtime_error("ggcat color runs do not match the number of k-mers in a unitig");
        }
        color_run_count += runs.size();

        ggcat_add_unitig_boundary_endpoints(endpoint_records,
                                            ks->seq.s,
                                            len,
                                            param.k,
                                            runs.front().first,
                                            runs.back().first);

        for (size_t r = 0; r < runs.size(); ++r) {
            uint32_t sid = runs[r].first;
            uint32_t count = runs[r].second;
            if (count > 1) {
                uint64_t c_count = color_cache.cardinality(sid);
                if (c_count > 0) hist[ggcat_hist_idx(c_count, c_count)] += (uint64_t)(count - 1);
                internal_edges += (uint64_t)(count - 1);
            }
            if (r + 1 < runs.size() && count > 0 && runs[r + 1].second > 0) {
                uint32_t next_sid = runs[r + 1].first;
                uint64_t c_count = 0;
                if (sid == next_sid) {
                    c_count = color_cache.cardinality(sid);
                } else {
                    c_count = intersect_size_sorted(color_cache.colors(sid), color_cache.colors(next_sid));
                    ++mixed_internal_edges;
                }
                if (c_count > 0) ++hist[ggcat_hist_idx(c_count, c_count)];
                ++internal_edges;
            }
        }
        ++unitigs;
    }
    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr, "Read %" PRIu64 " ggcat unitigs\n", unitigs);
    fprintf(stderr, "Color runs: %" PRIu64 "\n", color_run_count);
    fprintf(stderr, "Internal edge records: %" PRIu64 "\n", internal_edges);
    fprintf(stderr, "Mixed-subset internal edge records: %" PRIu64 "\n", mixed_internal_edges);
    fprintf(stderr, "Distinct color subsets: %zu\n", subset_ids.size());
    fprintf(stderr, "Endpoint records: %zu\n", endpoint_records.size());

    std::sort(endpoint_records.begin(), endpoint_records.end(), ggcat_endpoint_record_less);
    endpoint_records.erase(std::unique(endpoint_records.begin(), endpoint_records.end(), ggcat_endpoint_record_equal), endpoint_records.end());
    fprintf(stderr, "Unique endpoint records: %zu\n", endpoint_records.size());

    GgcatEndpointScratch scratch(n_colors);
    size_t groups = 0;
    size_t begin = 0;
    while (begin < endpoint_records.size()) {
        size_t end = begin + 1;
        while (end < endpoint_records.size() && endpoint_records[end].infix == endpoint_records[begin].infix) ++end;
        ggcat_add_endpoint_group_to_hist(endpoint_records, begin, end, param.k, color_cache, scratch, hist);
        ++groups;
        begin = end;
    }
    fprintf(stderr, "Endpoint infix groups: %zu\n", groups);

    for (size_t i = 0; i < hist.size(); ++i) printf("%" PRIu64 "\n", hist[i]);
#endif
}
