#ifndef TGUTILS_HPP
#define TGUTILS_HPP

#include <sstream>
#include <string>
#include <iterator>
#include <vector>
#include <cctype>

auto removeWhitespaces = [](std::string &line) -> void {
	line.erase( line.begin(), std::find_if_not( line.begin(), line.end(), [](int c) {return isspace(c);} ) ); // clean leading whitespaces
	line.erase( std::find_if_not( line.rbegin(),line.rend(), [](int c) {return isspace(c);} ).base(), line.end() ); // clean trailing whitespaces
};

class FileLine {
	std::string data;
	public:
		friend std::istream& operator>>(std::istream& is, FileLine& l) {
			return std::getline(is, l.data);
		}
		operator std::string() const {
			return data;
		}
};

template<typename T>
std::string JSONPair(const std::string &key, const std::vector<T> &value) {
	std::stringstream output; output << "\"" << key << "\":[";
	std::copy(value.begin(), value.end(), std::ostream_iterator<T>(output, ", "));
	std::string tmpout = output.str();

	unsigned int end = value.size() > 0 ? tmpout.size()-2 : tmpout.size();
	return tmpout.substr(0, end) + "]";
}

template<typename T>
std::string JSONPair(const std::string &key, const T &value) {
	std::stringstream output; output << "\"" << key << "\":" << value;
	return output.str();
}

std::string JSONPair(const std::string &key, const std::string &value);
std::string JSONPair(const std::string &key, const bool &value);
std::string JSONPair(const std::string &key, const char &value);
std::string JSON(const std::vector<std::string> &KVpairs);

std::string toString(char c);

double vectorDistance(const std::vector<double> &v1, const std::vector<double> &v2);
double squaredVectorDistance(const std::vector<double> &v1, const std::vector<double> &v2);

#endif
