#ifndef __PARSEHEADER_HH__
#define __PARSEHEADER_HH__

#include <stdio.h>
#include <string>
#include <filesystem>
#include <vector>
#include <unordered_map>
#include <variant>

#include "detail/phDriver.hh"

namespace fs = std::filesystem;

#define MUST_DEFINE true
#define DONT_CARE false

class HeaderStream {
public:
    HeaderStream(const fs::path &fn);
    HeaderStream(const char *in_buffer, size_t in_bufferlength, const fs::path &source_name);
    virtual ~HeaderStream(void);

    void OpenForRead(void);
    void Close(void);
    void SkipHeader(void);
    static size_t SkipHeaderFromFP(FILE *);
    void ReadHeader(void);

    fs::path name;
    char *buffer;
    size_t bufferlength;
    FILE *fp;

private:
    bool use_memory_buffer = false;
    bool owns_buffer = true;

    void GetHeaderLength(void);
};

void OpenStreamForWrite(std::ofstream& stream, const fs::path &fn, bool overwrite);
FILE* OpenForWrite(const fs::path &fn, bool overwrite);
void WriteHStream(FILE *fp, const std::string &m);
void WriteHStream(FILE *fp, const std::string &m, const std::string &pre);
void WriteHStream(FILE *fp, HeaderStream &in);
void WriteHStream(FILE *fp, HeaderStream &in, const std::string &pre);
void FinalizeHeader(FILE *fout);
void FinalizeHeader(std::stringstream &ss);

class phDriver;

class ParseHeader {
public:
    ParseHeader(void);
    ~ParseHeader() = default;

    // Disable copy and move, as usually this is not what the user wants.
    ParseHeader(const ParseHeader&) = delete;
    ParseHeader& operator=(const ParseHeader&) = delete;
    ParseHeader(ParseHeader&&) = delete;
    ParseHeader& operator=(ParseHeader&&) = delete;

    // register variables with the parser
    template <typename T>
    void register_vars(T &param);

    // install a scalar
    template <typename T>
    void installscalar(const std::string &name, T& var, bool must_define);

    // Install a vector
    template <typename T>
    void installvector(const std::string &name, std::vector<T> &var, bool must_define, size_t maxlen = 1024);

    void ReadHeader(HeaderStream &in);
    void ParseBuffer(HeaderStream &in);

private:
    std::unique_ptr<phDriver> phdriver;

    using VectorPtr = std::variant<
        std::vector<int> *,
        std::vector<double> *,
        std::vector<std::string> *,
        std::vector<fs::path> *
    >;

    struct VectorBinding {
        size_t original_size; // size before we over-allocate
        VectorPtr vec;
    };

    std::unordered_map<std::string, VectorBinding> vectors;

    void resize_vectors(void);
};

template<typename T>
inline void ParseHeader::installscalar(const std::string &name, T &var, bool must_define) {
    phdriver->InstallSym(name, &var, 1, 0, must_define);
}

template <typename T>
inline void ParseHeader::installvector(const std::string &name, std::vector<T> &var, bool must_define, size_t maxlen) {
    const size_t original_size = var.size();

    // Don't accidentally truncate user defaults if they already sized the vector.
    if (maxlen < original_size) maxlen = original_size;

    var.resize(maxlen);
    phdriver->InstallSym(name, var.data(), maxlen, 1, must_define);

    vectors[name] = VectorBinding{original_size, &var};  // will resize after parsing
}

#endif // __PARSEHEADER_HH__