#pragma once

#include <zlib.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
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

static inline bool ggcat_encode_window(const char *seq, int pos, int len, uint64_t *out) {
    uint64_t x = 0;
    for (int i = 0; i < len; ++i) {
        int c = seq_nt4_table[(uint8_t)seq[pos + i]];
        if (c >= 4) return false;
        x = (x << 2) | (uint64_t)c;
    }
    *out = x;
    return true;
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

struct GgcatEndpointGroup {
    std::array<std::vector<uint32_t>, 4> left;
    std::array<std::vector<uint32_t>, 4> right;
};

static inline void ggcat_push_unique(std::vector<uint32_t>& xs, uint32_t x) {
    if (std::find(xs.begin(), xs.end(), x) == xs.end()) xs.push_back(x);
}

static inline void ggcat_add_left_endpoint(std::unordered_map<uint64_t, GgcatEndpointGroup>& groups,
                                            uint64_t kmer, int k, uint32_t sid) {
    uint64_t infix = kmer & ggcat_mask_bits(k - 1);
    uint8_t a = (uint8_t)((kmer >> (2 * (k - 1))) & 3ULL);
    uint64_t rc = ggcat_revcomp_bits(infix, k - 1);
    if (infix < rc) {
        ggcat_push_unique(groups[infix].left[a], sid);
    } else {
        ggcat_push_unique(groups[rc].right[ggcat_comp_base(a)], sid);
    }
}

static inline void ggcat_add_right_endpoint(std::unordered_map<uint64_t, GgcatEndpointGroup>& groups,
                                             uint64_t kmer, int k, uint32_t sid) {
    uint64_t infix = kmer >> 2;
    uint8_t b = (uint8_t)(kmer & 3ULL);
    uint64_t rc = ggcat_revcomp_bits(infix, k - 1);
    if (infix < rc) {
        ggcat_push_unique(groups[infix].right[b], sid);
    } else {
        ggcat_push_unique(groups[rc].left[ggcat_comp_base(b)], sid);
    }
}

static inline uint64_t ggcat_pair_key(uint32_t a, uint32_t b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)a << 32) | (uint64_t)b;
}

static inline uint64_t ggcat_hist_idx(uint64_t j, uint64_t sigma) {
    return (sigma * (sigma - 1)) / 2 + j - 1;
}

static std::vector<std::pair<uint32_t, uint32_t> > ggcat_parse_color_runs(const std::string& header) {
    std::vector<std::pair<uint32_t, uint32_t> > runs;
    size_t pos = 0;
    while ((pos = header.find("C:", pos)) != std::string::npos) {
        pos += 2;
        char *end = 0;
        uint32_t sid = (uint32_t)strtoul(header.c_str() + pos, &end, 16);
        if (end == header.c_str() + pos || *end != ':') continue;
        uint32_t count = (uint32_t)strtoul(end + 1, &end, 10);
        runs.push_back(std::make_pair(sid, count));
        pos = (size_t)(end - header.c_str());
    }
    return runs;
}

static std::vector<uint32_t> ggcat_expand_runs(const std::vector<std::pair<uint32_t, uint32_t> >& runs, size_t n) {
    std::vector<uint32_t> sids;
    sids.reserve(n);
    for (size_t i = 0; i < runs.size(); ++i) {
        for (uint32_t j = 0; j < runs[i].second; ++j) sids.push_back(runs[i].first);
    }
    if (sids.size() != n) {
        throw std::runtime_error("ggcat color runs do not match the number of k-mers in a unitig");
    }
    return sids;
}

static inline void ggcat_union_intersection(const std::vector<uint32_t>& a,
                                            const std::vector<uint32_t>& b,
                                            std::vector<uint32_t>& out) {
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            out.push_back(a[i]);
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
}

static inline void ggcat_sort_unique(std::vector<uint32_t>& xs) {
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
}

static void ggcat_add_endpoint_group_to_hist(uint64_t infix,
                                             const GgcatEndpointGroup& group,
                                             int k,
                                             GgcatColorCache& colors,
                                             std::vector<uint64_t>& hist) {
    std::array<std::vector<uint32_t>, 16> edge_colors;

    for (uint8_t a = 0; a < 4; ++a) {
        for (uint8_t b = 0; b < 4; ++b) {
            uint64_t kp1mer = ((uint64_t)a << (2 * k)) | (infix << 2) | (uint64_t)b;
            GgcatCanonicalEdge ce;
            if (!ggcat_canonical_edge(kp1mer, k, &ce)) continue;
            if (ce.infix != infix || ce.edge != (uint8_t)((a << 2) | b)) continue;

            std::vector<uint32_t>& dst = edge_colors[(a << 2) | b];
            for (size_t li = 0; li < group.left[a].size(); ++li) {
                const std::vector<uint32_t>& lc = colors.colors(group.left[a][li]);
                for (size_t ri = 0; ri < group.right[b].size(); ++ri) {
                    const std::vector<uint32_t>& rc = colors.colors(group.right[b][ri]);
                    ggcat_union_intersection(lc, rc, dst);
                }
            }
            ggcat_sort_unique(dst);
        }
    }

    std::unordered_map<uint32_t, uint8_t> multiplicity;
    for (int e = 0; e < 16; ++e) {
        for (size_t i = 0; i < edge_colors[e].size(); ++i) {
            uint8_t& m = multiplicity[edge_colors[e][i]];
            if (m < 2) ++m;
        }
    }

    uint64_t sigma = multiplicity.size();
    if (sigma == 0) return;

    for (int e = 0; e < 16; ++e) {
        uint64_t j = 0;
        for (size_t i = 0; i < edge_colors[e].size(); ++i) {
            if (multiplicity[edge_colors[e][i]] == 1) ++j;
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

    std::vector<uint64_t> hist((n_colors + 1) * n_colors / 2, 0);
    std::unordered_map<uint64_t, uint64_t> internal_pairs;
    std::unordered_map<uint64_t, GgcatEndpointGroup> endpoint_groups;

    gzFile fp = gzopen(graph_path, "r");
    if (fp == 0) {
        fprintf(stderr, "Could not open ggcat graph: %s\n", graph_path);
        exit(EXIT_FAILURE);
    }
    kseq_t *ks = kseq_init(fp);
    int ret;
    uint64_t unitigs = 0, internal_edges = 0;
    while ((ret = kseq_read(ks)) >= 0) {
        int len = (int)ks->seq.l;
        if (len < param.k) continue;
        size_t n_kmers = (size_t)(len - param.k + 1);

        std::string header;
        if (ks->name.l) header.append(ks->name.s, ks->name.l);
        if (ks->comment.l) {
            header.push_back(' ');
            header.append(ks->comment.s, ks->comment.l);
        }
        std::vector<std::pair<uint32_t, uint32_t> > runs = ggcat_parse_color_runs(header);
        if (runs.empty()) {
            throw std::runtime_error("ggcat unitig has no color runs; build the graph with ggcat build -c");
        }
        std::vector<uint32_t> sids = ggcat_expand_runs(runs, n_kmers);

        uint64_t first_kmer = 0, last_kmer = 0;
        if (ggcat_encode_window(ks->seq.s, 0, param.k, &first_kmer)) {
            ggcat_add_right_endpoint(endpoint_groups, first_kmer, param.k, sids.front());
        }
        if (ggcat_encode_window(ks->seq.s, (int)n_kmers - 1, param.k, &last_kmer)) {
            ggcat_add_left_endpoint(endpoint_groups, last_kmer, param.k, sids.back());
        }

        for (size_t i = 0; i + 1 < n_kmers; ++i) {
            uint64_t kp1mer = 0;
            if (!ggcat_encode_window(ks->seq.s, (int)i, param.k + 1, &kp1mer)) continue;
            GgcatCanonicalEdge ce;
            if (!ggcat_canonical_edge(kp1mer, param.k, &ce)) continue;
            ++internal_pairs[ggcat_pair_key(sids[i], sids[i + 1])];
            ++internal_edges;
        }
        ++unitigs;
    }
    kseq_destroy(ks);
    gzclose(fp);

    fprintf(stderr, "Read %" PRIu64 " ggcat unitigs\n", unitigs);
    fprintf(stderr, "Internal edge records: %" PRIu64 "\n", internal_edges);
    fprintf(stderr, "Endpoint infix groups: %zu\n", endpoint_groups.size());

    for (std::unordered_map<uint64_t, uint64_t>::const_iterator it = internal_pairs.begin(); it != internal_pairs.end(); ++it) {
        uint32_t sid1 = (uint32_t)(it->first >> 32);
        uint32_t sid2 = (uint32_t)(it->first & 0xffffffffULL);
        uint64_t c_count = intersect_size_sorted(color_cache.colors(sid1), color_cache.colors(sid2));
        if (c_count > 0) hist[ggcat_hist_idx(c_count, c_count)] += it->second;
    }

    for (std::unordered_map<uint64_t, GgcatEndpointGroup>::const_iterator it = endpoint_groups.begin(); it != endpoint_groups.end(); ++it) {
        ggcat_add_endpoint_group_to_hist(it->first, it->second, param.k, color_cache, hist);
    }

    for (size_t i = 0; i < hist.size(); ++i) printf("%" PRIu64 "\n", hist[i]);
#endif
}
