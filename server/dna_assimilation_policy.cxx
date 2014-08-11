/*
 * Copyright 2012, 2009 Travis Desell and the University of North Dakota.
 *
 * This file is part of the Toolkit for Asynchronous Optimization (TAO).
 *
 * TAO is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TAO is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TAO.  If not, see <http://www.gnu.org/licenses/>.
 * */

#include <vector>
#include <cstdlib>
#include <string>
#include <fstream>

#include <math.h>

#include "config.h"
#include "util.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "md5_file.h"
#include "error_numbers.h"
#include "validate_util.h"

#include "stdint.h"
#include "mysql.h"
#include "boinc_db.h"

#include "undvc_common/file_io.hxx"
#include "undvc_common/parse_xml.hxx"

#include <boost/algorithm/string.hpp>

using namespace std;

#define mysql_query_check(conn, query) __mysql_check (conn, query, __FILE__, __LINE__)

void __mysql_check(MYSQL *conn, string query, const char *file, const int line) {
    mysql_query(conn, query.c_str());

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR in MySQL query: '" << query.c_str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << file << ":" << line;
        cerr << ex_msg.str() << endl;
        exit(1);
    }   
}

MYSQL *wildlife_db_conn = NULL;

void initialize_database() {
    wildlife_db_conn = mysql_init(NULL);

    //shoud get database info from a file
    string db_host, db_name, db_password, db_user;
    ifstream db_info_file("../wildlife_db_info");

    db_info_file >> db_host >> db_name >> db_user >> db_password;
    db_info_file.close();

    cout << "parsed db info:" << endl;
    cout << "\thost: " << db_host << endl;
    cout << "\tname: " << db_name << endl;
    cout << "\tuser: " << db_user << endl;
    cout << "\tpass: " << db_password << endl;

    if (mysql_real_connect(wildlife_db_conn, db_host.c_str(), db_user.c_str(), db_password.c_str(), db_name.c_str(), 0, NULL, 0) == NULL) {
        cerr << "Error connecting to database: " << mysql_errno(wildlife_db_conn) << ", '" << mysql_error(wildlife_db_conn) << "'" << endl;
        exit(1);
    }
}


//returns 0 on sucess
int assimilate_handler(WORKUNIT& wu, vector<RESULT>& results, RESULT& canonical_result) {
    if (wildlife_db_conn == NULL) initialize_database();

    //need to read wu.xml_doc
    string xml_doc;

    ostringstream oss;
    oss << "SELECT xml_doc FROM workunit WHERE id = " << wu.id;
    string query = oss.str();

    mysql_query(boinc_db.mysql, query.c_str());

    MYSQL_RES *my_result = mysql_store_result(boinc_db.mysql);
    if (mysql_errno(boinc_db.mysql) == 0) {
        MYSQL_ROW row = mysql_fetch_row(my_result);

        if (row == NULL) {
            log_messages.printf(MSG_CRITICAL, "Could not get row from workunit with query '%s'. Error: %d -- '%s'\n", xml_doc.c_str(), mysql_errno(boinc_db.mysql), mysql_error(boinc_db.mysql));
            exit(1);
        }

        xml_doc = row[0];
    } else {
        log_messages.printf(MSG_CRITICAL, "Could execute query '%s'. Error: %d -- '%s'\n", xml_doc.c_str(), mysql_errno(boinc_db.mysql), mysql_error(boinc_db.mysql));
        exit(1);
    }
    mysql_free_result(my_result);

    /*
     * Now that the workunit xml has been gotten, we can parse it for the appropriate information
     */
    string tag_str;
    
    try {
        tag_str = parse_xml<string>(xml_doc, "tag");
    } catch (string error_message) {
        log_messages.printf(MSG_CRITICAL, "wildlife_assimilation_policy assimilate_handler([RESULT#%d %s]) failed with error: %s\n", canonical_result.id, canonical_result.name, error_message.c_str());
        log_messages.printf(MSG_CRITICAL, "XML:\n'%s'\n", xml_doc.c_str());
        exit(1);

        return 0;
    }

    vector<double> p1;
    try {
        string prob_str = parse_xml<string>(canonical_result.stderr_out, "slice_probabilities");

        istringstream iss1(prob_str);
        copy(istream_iterator<double>(iss1), istream_iterator<double>(), back_inserter<vector<double> >(p1));

    } catch (string error_message) {
        log_messages.printf(MSG_CRITICAL, "wildlife_assimilation_policy assimilate_handler([RESULT#%d %s]) failed with error: %s\n", canonical_result.id, canonical_result.name, error_message.c_str());
        log_messages.printf(MSG_CRITICAL, "XML:\n%s\n", canonical_result.stderr_out);

        return 0;
    }

    cout << "tag: " << tag_str << endl;
    cout << "result name: " << canonical_result.name << endl;

    cout << "parsed probabilities: " << p1.size() << endl;
    for (uint32_t i = 0; i < p1.size(); i++) {
        cout << "\t" << p1[i] << endl;
    }

    //get video id
    //result name is video_<video_id>_<time>_<result number>
    string result_name = results[0].name;
/*    string result_name = canonical_result.name;     SHOULD BE THIS */
    uint32_t first_pos = result_name.find("_", 0) + 1;
    uint32_t second_pos = result_name.find("_", first_pos);

    if (first_pos == string::npos || second_pos == string::npos) {
        log_messages.printf(MSG_CRITICAL, "wildlife_assimilation_policy assimilate_handler failed with 'malformed result name error', result name: %s\n", result_name.c_str());
        exit(1);
    }

    string video_id = result_name.substr( first_pos, (second_pos - first_pos) );

    cout << "parsed video id: '" << video_id << "'" << endl;

    //get the video segment id, species_id and location_id:
    //  SELECT id, species_id, location_id FROM video_segment_2 WHERE video_id = video_id and number = i

    ostringstream full_video_query;
    full_video_query << "SELECT species_id, location_id, start_time, duration_s FROM video_2 WHERE id = " << video_id << endl;

    mysql_query_check(wildlife_db_conn, full_video_query.str());
    MYSQL_RES *video_result = mysql_store_result(wildlife_db_conn);

    cout << " got video result" << endl;

    MYSQL_ROW full_video_row = mysql_fetch_row(video_result);
    int species_id = atoi(full_video_row[0]);
    int location_id = atoi(full_video_row[1]);
    string start_time = full_video_row[2];
    int duration_s = atoi(full_video_row[3]);

    uint32_t calculated_segments = ceil(duration_s / 180.0);
    cout << "duration_s: " << duration_s << ", segments: " << calculated_segments << ", p1.size(): " << p1.size() << endl;

    if (p1.size() < (calculated_segments - 1) || p1.size() > calculated_segments) {
        cout << "Wrong number of seconds" << endl;
        mysql_free_result(video_result);
        exit(1);
        return 0;
    }

    for (uint32_t i = 0; i < p1.size(); i++) {
        //determine video segment id for each 3 minute probability
        ostringstream video_segment_query;
        video_segment_query << "SELECT id FROM video_segment_2 WHERE video_id = " << video_id << " AND number = " << i;

        cout << "getting video segments: '" << video_segment_query.str() << "'" << endl;

        mysql_query_check(wildlife_db_conn, video_segment_query.str());
        MYSQL_RES *video_segment_result = mysql_store_result(wildlife_db_conn);

        if (video_segment_result == NULL) {
            log_messages.printf(MSG_CRITICAL, "wildlife_assimilation_policy assimilate_handler failed with 'no matching video segment id', result name: %s\n", result_name.c_str());
            log_messages.printf(MSG_CRITICAL, "\tvideo segment id = ? for video_id = %s and number = %u\n", video_id.c_str(), i);
            exit(1);
        }

        MYSQL_ROW video_segment_row = mysql_fetch_row(video_segment_result);

        if (!video_segment_row) {
            cout << "could not get a video segment row" << endl;
            exit(1);
        }

        int video_segment_id = atoi(video_segment_row[0]);

        cout << "video segment id = " << video_segment_id << " for video_id = " << video_id << " and number = " << i << endl;


        //Set the data for the classification for that three minute segment
        ostringstream classifications_query;
        classifications_query << "REPLACE INTO classifications "
                              << " SET video_id = " << video_id
                              << ", video_segment_id = " << video_segment_id
                              << ", probability = " << p1[i]
//                              << ", type = '" << detection_method << "'"
//                              << ", detection = '" << detection_target << "'"
                              << ", tag = '" << tag_str << "'"
                              << ", species_id = " << species_id
                              << ", location_id = " << location_id
                              << ", start_time = (SELECT ADDTIME(\"" << start_time << "\", SEC_TO_TIME(" << i * 180 << ")))";

        //this will allow multiple classifications with the same information

        cout << "classifications query: '" << classifications_query.str() << "'" << endl;

        mysql_query_check(wildlife_db_conn, classifications_query.str());
        MYSQL_RES *classifications_result = mysql_store_result(wildlife_db_conn);

        /*
        if (classifications_result == NULL) {
            log_messages.printf(MSG_CRITICAL, "wildlife_assimilation_policy assimilate_handler failed with 'could not insert classification', result name: %s\n", result_name.c_str());
            exit(1);
        }
        */

        cout << "freeing results" << endl;

        mysql_free_result(video_segment_result);
        mysql_free_result(classifications_result);

        cout << "next iteration" << endl;
    }

    cout << "finished loop" << endl;

    mysql_free_result(video_result);

    //Don't need to do anything, when the result is validated it gets inserted into the database directly
    return 0;
}
