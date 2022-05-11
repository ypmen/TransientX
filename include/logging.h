/**
 * @author Yunpeng Men
 * @email ypmen@mpifr-bonn.mpg.de
 * @create date 2021-06-16 10:24:04
 * @modify date 2021-06-16 10:24:04
 * @desc [description]
 */

#ifndef LOGGING_H
#define LOGGING_H

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <sstream>
#include <iomanip>

namespace logging = boost::log;
inline void init_logging()
{
	logging::core::get()->set_filter
	(
		logging::trivial::severity >= logging::trivial::info
	);
}

inline void format_logging(const std::string &title, const std::vector<std::pair<std::string, std::string>> &meta)
{
	std::stringstream s_meta;
	for (auto it=meta.begin(); it!=meta.end(); ++it)
	{
		s_meta<<std::left<<std::setw(32)<<it->first<<":"<<std::setw(31)<<it->second<<'\n';
	}

	BOOST_LOG_TRIVIAL(info) << '\n'
	<< std::setfill('=') << std::setw((64+title.size())/2) << title << std::setw((64-title.size())/2) << "" << '\n'
	<< s_meta.str()
	<< std::setfill('=') << std::setw(64) << "";
}

#endif /* LOGGING_H */
