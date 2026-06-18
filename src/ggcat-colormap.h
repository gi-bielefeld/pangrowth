#pragma once

#include <stdint.h>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef PANGROWTH_WITH_GGCAT
#include <ggcat.hh>
static void pangrowth_ggcat_silent_messages(ggcat::MessageLevel, const char*) {}
#endif

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
        num_colors_ = ggcat::GGCATInstance::dump_colors(colormap_path_).size();
#else
        (void)threads_;
        num_colors_ = 0;
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
        return colors(subset_id).size();
    }

    size_t num_colors() const {
        return num_colors_;
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
                std::vector<uint32_t> v(cols.data, cols.data + cols.size);
                std::sort(v.begin(), v.end());
                v.erase(std::unique(v.begin(), v.end()), v.end());
                cache_[subset].swap(v);
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
    std::unordered_map<uint32_t, std::vector<uint32_t> > cache_;
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
