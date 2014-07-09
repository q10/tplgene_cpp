#include <cmath>
#include "TGUtils.hpp"

using namespace std;

string JSONPair(const string &key, const string &value) {
	stringstream output; output << "\"" << key << "\":\"" << value << "\"";
	return output.str();
}

string JSONPair(const string &key, const bool &value) {
	stringstream output; output << "\"" << key << "\":" << (value ? "true" : "false");
	return output.str();
}

string JSONPair(const string &key, const char &value) {
	stringstream output; output << "\"" << key << "\":\'" << value << "\'";
	return output.str();
}

string JSON(const vector<string> &KVpairs) {
	stringstream output; output << "{";
	std::copy(KVpairs.begin(), KVpairs.end(), ostream_iterator<string>(output, ", "));
	string tmpout = output.str();

	unsigned int end = KVpairs.size() > 0 ? tmpout.size()-2 : tmpout.size();
	return tmpout.substr(0, end) + "}";
}

string toString(char c) {
	stringstream ss; ss << c; return ss.str();
}

double vectorDistance(const vector<double> &v1, const vector<double> &v2) {
	return sqrt( squaredVectorDistance(v1, v2) );
}

double squaredVectorDistance(const vector<double> &v1, const vector<double> &v2) {
	unsigned int size = min(v1.size(), v2.size());
	double sum = 0;
	for (unsigned int i=0; i < size; i++) {
		double tmp = v1[i] - v2[i];
		sum += tmp*tmp;
	} return sum;
}