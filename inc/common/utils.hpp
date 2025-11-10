#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <sstream>
#include <vector>

std::vector<std::string> split_string(const std::string& str, const std::string& delim);
std::string join(const std::vector<std::string>& v, const std::string& delim);
void print_info(const char* msg);
void print_warning(const char* msg);
void print_error(const char* msg);
void print_debug(const char* msg);
void parse_err_msg(int line, const char* msg);

bool file_exists(const std::string& name);
std::string read_file(const std::string& filename);
std::string read_until(std::ifstream& file, const std::string& delim, bool include_delim = true, std::string ignore_chars = "");

uint32_t log_2(const uint32_t x);

template<std::unsigned_integral T>
std::string uitoh(T x) {
    std::stringstream ss;
    ss << std::hex << x;
    return ss.str();
}
#endif