#pragma once

#include <stdint.h>
#include <cassert>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef PANGROWTH_WITH_GGCAT
#include <ggcat.hh>
static void pangrowth_ggcat_silent_messages(ggcat::MessageLevel, const char*) {}

#ifndef NDEBUG
static inline bool ggcat_colors_are_sorted_unique(ggcat::Slice<uint32_t> cols) {
    for (size_t i = 1; i < cols.size; ++i) {
        if (cols.data[i - 1] >= cols.data[i]) return false;
    }
    return true;
}
#endif
#endif

struct GgcatColormapHeader {
    size_t colors;
    size_t color_subsets;
};

static inline uint64_t ggcat_read_le64(const uint8_t *p) {
    uint64_t x = 0;
    for (int i = 7; i >= 0; --i) x = (x << 8) | (uint64_t)p[i];
    return x;
}

static inline GgcatColormapHeader ggcat_read_colormap_header(const std::string& colormap_path) {
    uint8_t header[56];
    FILE *fp = std::fopen(colormap_path.c_str(), "rb");
    if (fp == 0) throw std::runtime_error("could not open ggcat colormap file");
    size_t n = std::fread(header, 1, sizeof(header), fp);
    std::fclose(fp);
    if (n != sizeof(header)) throw std::runtime_error("could not read ggcat colormap header");
    if (std::memcmp(header, "GGCAT_CMAP_RNLEN", 16) != 0) {
        throw std::runtime_error("unsupported ggcat colormap format");
    }
    GgcatColormapHeader h;
    h.colors = (size_t)ggcat_read_le64(header + 32);
    h.color_subsets = (size_t)ggcat_read_le64(header + 40);
    return h;
}

class GgcatColorCache {
public:
    explicit GgcatColorCache(const std::string& graph_path, size_t threads = 1)
        : graph_path_(graph_path), threads_(threads) {
#ifdef PANGROWTH_WITH_GGCAT
        ggcat::GGCATConfig config;
        config.use_temp_dir = true;
        config.temp_dir = ".ggcat_tmp";
        config.memory = 2.0;
        config.prefer_memory = false;
        config.total_threads_count = threads_;
        config.intermediate_compression_level = static_cast<uint32_t>(-1);
        config.use_stats_file = false;
        config.stats_file = "";
        config.messages_callback = pangrowth_ggcat_silent_messages;
        instance_ = ggcat::GGCATInstance::create(config);
        colormap_path_ = ggcat::GGCATInstance::get_colormap_file(graph_path_);
        GgcatColormapHeader header = ggcat_read_colormap_header(colormap_path_);
        num_colors_ = header.colors;
        num_color_subsets_ = header.color_subsets;
#else
        (void)threads_;
        num_colors_ = 0;
        num_color_subsets_ = 0;
        throw std::runtime_error("pangrowth was built without ggcat API support");
#endif
    }

    const std::vector<uint32_t>& colors(uint32_t subset_id) {
        std::unordered_map<uint32_t, std::vector<uint32_t> >::iterator it = cache_.find(subset_id);
        if (it != cache_.end()) return it->second;
        query(&subset_id, 1);
        it = cache_.find(subset_id);
        if (it == cache_.end()) throw std::runtime_error("ggcat colormap query did not return requested subset");
        return it->second;
    }

    size_t cardinality(uint32_t subset_id) {
        std::unordered_map<uint32_t, std::vector<uint32_t> >::iterator it = cache_.find(subset_id);
        if (it != cache_.end()) return it->second.size();
        std::unordered_map<uint32_t, size_t>::iterator cit = cardinality_cache_.find(subset_id);
        if (cit != cardinality_cache_.end()) return cit->second;
        query_cardinalities(&subset_id, 1);
        cit = cardinality_cache_.find(subset_id);
        if (cit == cardinality_cache_.end()) throw std::runtime_error("ggcat colormap query did not return requested subset");
        return cit->second;
    }

    size_t num_colors() const {
        return num_colors_;
    }

    size_t num_color_subsets() const {
        return num_color_subsets_;
    }

    void reserve(size_t n) {
        cache_.reserve(n);
        cardinality_cache_.reserve(n);
    }

    void query(const uint32_t* subset_ids, size_t n) {
#ifdef PANGROWTH_WITH_GGCAT
        std::vector<uint32_t> missing;
        missing.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (cache_.find(subset_ids[i]) == cache_.end()) missing.push_back(subset_ids[i]);
        }
        if (missing.empty()) return;

        instance_->query_colormap(
            colormap_path_,
            missing.data(),
            missing.size(),
            true,
            [&](uint32_t subset, ggcat::Slice<uint32_t> cols) {
#ifndef NDEBUG
                assert(ggcat_colors_are_sorted_unique(cols));
#endif
                std::vector<uint32_t> v(cols.data, cols.data + cols.size);
                cache_[subset].swap(v);
            });
#else
        (void)subset_ids;
        (void)n;
        throw std::runtime_error("pangrowth was built without ggcat API support");
#endif
    }

    void query_cardinalities(const uint32_t* subset_ids, size_t n) {
#ifdef PANGROWTH_WITH_GGCAT
        std::vector<uint32_t> missing;
        missing.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            if (cache_.find(subset_ids[i]) == cache_.end()
                && cardinality_cache_.find(subset_ids[i]) == cardinality_cache_.end()) {
                missing.push_back(subset_ids[i]);
            }
        }
        if (missing.empty()) return;

        instance_->query_colormap(
            colormap_path_,
            missing.data(),
            missing.size(),
            true,
            [&](uint32_t subset, ggcat::Slice<uint32_t> cols) {
                cardinality_cache_[subset] = cols.size;
            });
#else
        (void)subset_ids;
        (void)n;
        throw std::runtime_error("pangrowth was built without ggcat API support");
#endif
    }

private:
    std::string graph_path_;
    std::string colormap_path_;
    size_t threads_;
    size_t num_colors_;
    size_t num_color_subsets_;
    std::unordered_map<uint32_t, std::vector<uint32_t> > cache_;
    std::unordered_map<uint32_t, size_t> cardinality_cache_;
#ifdef PANGROWTH_WITH_GGCAT
    ggcat::GGCATInstance* instance_;
#endif
};

static inline size_t intersect_size_sorted(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
    size_t i = 0, j = 0, n = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) {
            ++n;
            ++i;
            ++j;
        } else if (a[i] < b[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    return n;
}
