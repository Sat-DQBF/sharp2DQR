#include "common/utils.hpp"

#include <assert.h>
#include <fstream>
#include <sys/stat.h>
#include <vector>

// https://stackoverflow.com/questions/289347/using-strtok-with-a-stdstring
std::vector<std::string> split_string(const std::string& str, const std::string& delim) {
    std::vector<std::string> parts;
    size_t start, end = 0;
    while (end < str.size()) {
        start = end;
        while (start < str.size() && (delim.find(str[start]) != std::string::npos)) {
            start++;  // skip initial whitespace
        }
        end = start;
        while (end < str.size() && (delim.find(str[end]) == std::string::npos)) {
            end++;  // skip to end of word
        }
        if (end - start != 0) {  // just ignore zero-length strings.
            parts.push_back(std::string(str, start, end - start));
        }
    }
    return parts;
}

std::string join(const std::vector<std::string>& v, const std::string& delim) {
    std::string s;
    for (std::vector<std::string>::const_iterator it = v.begin(); it != v.end(); it++) {
        s += *it;
        if (it != v.end() - 1) {
            s += delim;
        }
    }
    return s;
}

void print_info(const char* msg) {
    std::vector<std::string> parts;
    parts = split_string(msg, "\n");
    for (auto it = parts.begin(); it < parts.end(); it++) {
        printf("\033[96m[INFO]\033[0m %s\n", (*it).c_str());
    }
    fflush(stdout);
}

void print_warning(const char* msg) {
    std::vector<std::string> parts;
    parts = split_string(msg, "\n");
    for (auto it = parts.begin(); it < parts.end(); it++) {
        printf("\033[93m[WARN]\033[0m %s\n", (*it).c_str());
    }
    fflush(stdout);
}

void print_error(const char* msg) {
    std::vector<std::string> parts;
    parts = split_string(msg, "\n");
    for (auto it = parts.begin(); it < parts.end(); it++) {
        printf("\033[91m[ERROR]\033[0m %s\n", (*it).c_str());
    }
    exit(-1);
}

void print_debug(const char *msg) {
    #ifdef DEBUG
    std::vector<std::string> parts;
    parts = split_string(msg, "\n");
    for (auto it = parts.begin(); it < parts.end(); it++) {
        printf("\033[92m[DEBUG]\033[0m %s\n", (*it).c_str());
    }
    fflush(stdout);
    #endif
}

void parse_err_msg(int line, const char* msg) {
    printf("\033[91m[ERROR]\033[0m Error parsing file at line %u: %s\n", line, msg);
    exit(-1);
}

bool file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

std::string read_file(const std::string &filename) {
    if (!file_exists(filename)) {
        print_error(("File " + filename + " does not exist").c_str());
    }
    std::ifstream file(filename);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

std::string read_until(std::ifstream &file, const std::string &delim, bool include_delim, std::string ignore_chars) {
    std::string s;
    assert (delim.find_first_of(ignore_chars) == std::string::npos);
    char c;
    int match_pos = 0;
    while (file.get(c)) {
        if (ignore_chars.find(c) != std::string::npos) {
            continue;
        }
        s += c;
        if (c == delim[match_pos]) {
            if (++match_pos == delim.size()) {
                break;
            }
        } else {
            match_pos = 0;
        }
    }
    if (!include_delim) {
        return s.substr(0, s.size() - delim.size());
    }
    return s;
}

uint32_t log_2(const uint32_t x) {
    return x == 0 ? 0 : 31 - __builtin_clz(x);
}