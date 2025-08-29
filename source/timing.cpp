/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugin timing
 *
 ******************************************************************************/

#include "timing.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <filesystem>
#include <system_error>
#include <limits>

using namespace std;
namespace Manta {

TimingData::TimingData() : updated(false), num(0) {
}

void TimingData::start(FluidSolver* parent, const string& name) {
	mLastPlugin = name;
	mPluginTimer.get();
}

void TimingData::stop(FluidSolver* parent, const string& name) {
	if (mLastPlugin == name && name != "FluidSolver::step") {
		updated = true;
		const string parentName = parent ? parent->getName() : "";
		MuTime diff = mPluginTimer.update();
		vector<TimingSet>& cur = mData[name];
		for (vector<TimingSet>::iterator it = cur.begin(); it != cur.end(); it++) {
			if (it->solver == parentName) {
				it->cur += diff;
				it->updated = true;
				return;
			}
		}
		TimingSet s;
		s.solver = parentName;
		s.cur = diff;
		s.updated = true;
		cur.push_back(s);
	}
}

void TimingData::step() {
	if (updated)
		num++;
	std::map<std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++) {
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (it2->updated) {
				it2->total += it2->cur;
				it2->num++;
			}
			it2->cur.clear();
			it2->updated = false;
		}
	}
	updated = false;
}
 
void TimingData::print() {
	MuTime total;
	total.clear();
	std::map<std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++)
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
			total += it2->cur;

	printf("\n-- STEP %3d ----------------------------\n", num);
	for (it = mData.begin(); it != mData.end(); it++) {
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (!it2->updated) continue;
			string name = it->first;
			if (it->second.size() > 1 && !it2->solver.empty())
				name += "[" + it2->solver + "]";
			printf("[%4.1f%%] %s (%s)\n", 100.0*((Real)it2->cur.time / (Real)total.time),
										  name.c_str(), it2->cur.toString().c_str());
		}
	}
	step();
		
	printf("----------------------------------------\n");
	printf("Total : %s\n\n", total.toString().c_str());
}

void TimingData::saveMean(const string& filename) {
	ofstream ofs(filename.c_str());
	step();
	if (!ofs.good())
		errMsg("can't open " + filename + " as timing log");
	ofs << "Mean timings of " << num << " steps :" <<endl <<endl;
	MuTime total;
	total.clear();
	std::map<std::string, std::vector<TimingSet> >::iterator it;
	for (it = mData.begin(); it != mData.end(); it++)
		for (vector<TimingSet>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			total += it2->cur;
			string name = it->first;
			if (it->second.size() > 1)
				name += "[" + it2->solver + "]";
			
			ofs << name << " " << (it2->total / it2->num) << endl;
		}
	 
	ofs << endl << "Total : " << total << " (mean " << total/num << ")" << endl;
	ofs.close();
}

static std::string escapeJsonString(const std::string &s) {
    std::ostringstream o;
    for (unsigned char c : s) {
        switch (c) {
            case '\"': o << "\\\""; break;
            case '\\': o << "\\\\"; break;
            case '\b': o << "\\b";  break;
            case '\f': o << "\\f";  break;
            case '\n': o << "\\n";  break;
            case '\r': o << "\\r";  break;
            case '\t': o << "\\t";  break;
            default:
                if (c < 0x20) {
                    o << "\\u"
                      << std::hex << std::setw(4) << std::setfill('0') << (int)c
                      << std::dec;
                } else {
                    o << c;
                }
        }
    }
    return o.str();
}

void TimingData::saveJson(const std::string& filename) {
    // compute total (same as print)
    MuTime total;
    total.clear();
    for (auto it = mData.begin(); it != mData.end(); ++it) {
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            total += it2->cur;
        }
    }

    std::ostringstream oss;
    oss << std::fixed; // we'll use fixed for percent formatting only when needed

    oss << "{\n";
    oss << "  \"step\": " << num << ",\n";
    oss << "  \"timings\": [\n";

    bool firstEntry = true;
    for (auto it = mData.begin(); it != mData.end(); ++it) {
        const std::string &groupName = it->first;
        for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
            if (!it2->updated) continue;

            if (!firstEntry) oss << ",\n";
            firstEntry = false;

            unsigned long ms = it2->cur.time; // keep as integer ms
            long double percent = 0.0L;
            if (total.time > 0) {
                percent = 100.0L * static_cast<long double>(ms) / static_cast<long double>(total.time);
            }

            oss << "    {\n";
            oss << "      \"name\": \"" << escapeJsonString(groupName) << "\",\n";
            if (!it2->solver.empty())
                oss << "      \"solver\": \"" << escapeJsonString(it2->solver) << "\",\n";
            oss << "      \"ms\": " << ms << ",\n";
            // percent with reasonable precision
            oss << "      \"percent\": " << std::setprecision(6) << static_cast<double>(percent) << "\n";
            oss << "    }";
            // stream state (precision/fixed) won't harm integer output later
        }
    }

    oss << "\n  ],\n";
    oss << "  \"total_ms\": " << total.time << "\n";
    oss << "}\n";

    // Ensure parent directory exists (C++17)
    std::error_code ec;
    std::filesystem::path p(filename);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path(), ec);
        if (ec) {
            std::cout << "Failed to create directory '" << p.parent_path().string()
                      << "': " << ec.message() << std::endl;
            return;
        }
    }

    // Write atomically to a temp file then rename
    std::string tmpname = p.string() + ".tmp";
    {
        std::ofstream ofs(tmpname, std::ios::binary);
        if (!ofs) {
            std::cout << "Error opening temp file for write: " << tmpname << std::endl;
            return;
        }
        ofs << oss.str();
        if (!ofs) {
            std::cout << "Error writing to temp file: " << tmpname << std::endl;
            return;
        }
        ofs.close();
    }

    std::filesystem::rename(tmpname, p, ec);
    if (ec) {
        std::cout << "Rename failed (" << ec.message() << "), attempting direct overwrite of '"
                  << p.string() << "'" << std::endl;
        // attempt fallback: remove temp and write directly
        std::remove(tmpname.c_str());
        std::ofstream ofs(p.string(), std::ios::binary);
        if (!ofs) {
            std::cout << "Error opening target file for write: " << p.string() << std::endl;
            return;
        }
        ofs << oss.str();
        if (!ofs) {
            std::cout << "Error writing to target file: " << p.string() << std::endl;
            return;
        }
        ofs.close();
    }

    // keep same behaviour as print(): advance step counter
    // step();
}
 
}

