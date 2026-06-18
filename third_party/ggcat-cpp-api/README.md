# ggcat C++ API wrapper

This directory vendors the small C++ wrapper from ggcat v2.0.0:

- `include/ggcat.hh`
- `include/ggcat-cpp-bindings.hh`
- `src/ggcat.cc`

The Rust ggcat static libraries are not vendored. Build them with:

```bash
./scripts/build_ggcat_cpp_api.sh
```

Then configure pangrowth with:

```bash
cmake -DPANGROWTH_WITH_GGCAT=ON \
  -DGGCAT_CPP_BINDINGS_LIB=$PWD/third_party/ggcat-cpp-api/lib/libggcat_cpp_bindings.a \
  -DGGCAT_CXX_INTEROP_LIB=$PWD/third_party/ggcat-cpp-api/lib/libggcat_cxx_interop.a \
  ..
```

The wrapper is MIT licensed; see `LICENSE`.
