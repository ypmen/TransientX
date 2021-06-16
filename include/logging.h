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

namespace logging = boost::log;
inline void init_logging()
{
    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::info
    );
}

#endif /* LOGGING_H */
