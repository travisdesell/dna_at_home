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
#include <unordered_map>

#include <math.h>
#include <sys/stat.h>
#include <sys/param.h>

#include "backend_lib.h"
#include "config.h"
#include "util.h"
#include "sched_config.h"
#include "sched_util.h"
#include "sched_msgs.h"
#include "md5_file.h"
#include "error_numbers.h"
#include "validate_util.h"
#include "str_replace.h"
#include "str_util.h"


#include "stdint.h"
#include "mysql.h"
#include "boinc_db.h"

//from undvc_common
#include "file_io.hxx"
#include "parse_xml.hxx"

#include <boost/algorithm/string.hpp>

using namespace std;

#define REPLICATION_FACTOR  2

#define mysql_query_check(query) __mysql_check (query, __FILE__, __LINE__)

//can just use boinc_db.mysql
void __mysql_check(string query, const char *file, const int line) {
    MYSQL *conn = boinc_db.mysql;
    mysql_query(conn, query.c_str());

    if (mysql_errno(conn) != 0) {
        ostringstream ex_msg;
        ex_msg << "ERROR in MySQL query: '" << query.c_str() << "'. Error: " << mysql_errno(conn) << " -- '" << mysql_error(conn) << "'. Thrown on " << file << ":" << line;
        cerr << ex_msg.str() << endl;
        exit(1);
    }   
}

bool file_exists(const string &filename) {
    struct stat buffer;   
    return (stat(filename.c_str(), &buffer) == 0); 
}

class GibbsSampler {
    public:
        int id;
        string name;
        string initial_workunit_filename;
        string checkpoint_workunit_filename;
        string result_xml_filename;
        string sequences_filename;
        int number_walks;
        int burn_in;
        int samples;
        int workunit_steps;
        string command_line_options;
        double rsc_fpops_bound;
        double rsc_fpops_est;
        double rsc_disk_bound;
        double rsc_memory_bound;
        int delay_bound;

        GibbsSampler(int id) {
            this->id = id;

            ostringstream oss;
            oss << "SELECT name, initial_workunit_xml_filename, checkpoint_workunit_xml_filename, result_xml_filename, sequences_filename, number_walks, burn_in, samples, workunit_steps, command_line_options, rsc_fpops_bound, rsc_fpops_est, rsc_disk_bound, rsc_memory_bound, delay_bound FROM gibbs_sampler WHERE id = " << id;
            mysql_query_check(oss.str());

            MYSQL_RES *sampler_result = mysql_store_result(boinc_db.mysql);
            MYSQL_ROW sampler_row = mysql_fetch_row(sampler_result);

            if (sampler_row == NULL) {
                log_messages.printf(MSG_CRITICAL, "ERROR: could not find gibbs sampler with id: %d\n", id);
                exit(1);
            }

            name = string(sampler_row[0]);
            initial_workunit_filename = string(sampler_row[1]);
            checkpoint_workunit_filename = string(sampler_row[2]);
            result_xml_filename = string(sampler_row[3]);
            sequences_filename = string(sampler_row[4]);
            number_walks = atoi(sampler_row[5]);
            burn_in = atoi(sampler_row[6]);
            samples = atoi(sampler_row[7]);
            workunit_steps = atoi(sampler_row[8]);
            command_line_options = string(sampler_row[9]);
            rsc_fpops_bound = atof(sampler_row[10]);
            rsc_fpops_est = atof(sampler_row[11]);
            rsc_disk_bound = atof(sampler_row[12]);
            rsc_memory_bound = atof(sampler_row[13]);
            delay_bound = atoi(sampler_row[14]);
        }

        string to_string() {
            ostringstream oss;

            oss << "[GIBBS SAMPLER"
                << " , id: " << id 
                << " , name: '" << name 
                << "', initial_workunit_filename = '" << initial_workunit_filename
                << "', checkpoint_workunit_filename = '" << checkpoint_workunit_filename
                << "', result_xml_filename = '" << result_xml_filename
                << "', seuences_filename = '" << sequences_filename
                << "', number_walks = " << number_walks
                << " , burn_in = " << burn_in
                << " , samples = " << samples
                << " , workunit_steps = " << workunit_steps
                << " , command_line_options = '" << command_line_options
                << "', rsc_fpops_bound = " << rsc_fpops_bound
                << " , rsc_fpops_est = " << rsc_fpops_est
                << " , rsc_disk_bound = " << rsc_disk_bound
                << " , rsc_memory_bound = " << rsc_memory_bound
                << " , delay_bound = " << delay_bound << "]";

            return oss.str();
        }
};


class GibbsWalk {
    public:
        int id;
        int sampler_id;
        int current_steps;
        string current_sites;
        int seed;

        GibbsWalk(int id) {
            this->id = id;

            ostringstream oss;
            oss << "SELECT sampler_id, current_steps, current_sites, seed FROM gibbs_walk WHERE id = " << id;
            mysql_query_check(oss.str());

            MYSQL_RES *walk_result = mysql_store_result(boinc_db.mysql);
            MYSQL_ROW walk_row = mysql_fetch_row(walk_result);

            if (walk_row == NULL) {
                log_messages.printf(MSG_CRITICAL, "ERROR: could not find gibbs walk with id: %d\n", id);
                exit(1);
            }
            sampler_id = atoi(walk_row[0]);
            current_steps = atoi(walk_row[1]);
            current_sites = string(walk_row[2]);
            seed = atoi(walk_row[3]);
        }

        string to_string() {
            ostringstream oss;

            oss << "[GIBBS WALK"
                << " , id: " << id 
                << " , sampler_id: '" << sampler_id
                << " , current_steps = " << current_steps
                << " , current_sites = '" << current_sites
                << "', seed = " << seed << "]";

            return oss.str();
        }

};

bool in_template_initialized = false;
char *in_template;

unordered_map<int, GibbsSampler*> samplers;

GibbsSampler* get_sampler(int id) {
    GibbsSampler *sampler;

    if (samplers.find(id) == samplers.end()) {
        sampler = new GibbsSampler(id);
        samplers[id] = sampler;
    } else {
        sampler = samplers[id];
    }

    return sampler;
}

//returns 0 on success
int assimilate_handler(WORKUNIT& wu, vector<RESULT>& results, RESULT& canonical_result) {
    //need to read wu.xml_doc
    ostringstream oss;
    oss << "SELECT xml_doc FROM workunit WHERE id = " << wu.id;
    string query = oss.str();

    mysql_query_check(query.c_str());

    MYSQL_RES *wu_result = mysql_store_result(boinc_db.mysql);
    MYSQL_ROW wu_row = mysql_fetch_row(wu_result);

    if (wu_row == NULL) {
        log_messages.printf(MSG_CRITICAL, "Could not get row from workunit with query '%s'. Error: %d -- '%s'\n", query.c_str(), mysql_errno(boinc_db.mysql), mysql_error(boinc_db.mysql));
        exit(1);
    }   

    string xml_doc = wu_row[0];

    mysql_free_result(wu_result);

    int sampler_id = -1;
    int walk_id = -1;
    try {
        sampler_id = parse_xml<int>(xml_doc, "sampler_id");
        walk_id = parse_xml<int>(xml_doc, "walk_id");
    } catch (string error_message) {
        log_messages.printf(MSG_CRITICAL, "dna_assimilator assimilate_handler([RESULT#%d %s]) failed with error: %s\n", canonical_result.id, canonical_result.name, error_message.c_str());
        log_messages.printf(MSG_CRITICAL, "XML:\n'%s'\n", xml_doc.c_str());
        exit(1);

        return 0;
    }   

    if (walk_id <= 0) {
        log_messages.printf(MSG_CRITICAL, "dna_assimilator assimilate_handler([RESULT#%d %s]) workunit had invalid walk_id (%d), should be > 0.\n", canonical_result.id, canonical_result.name, walk_id);
        return 0;
    }

    if (sampler_id <= 0) {
        log_messages.printf(MSG_CRITICAL, "dna_assimilator assimilate_handler([RESULT#%d %s]) workunit had invalid sampler_id (%d), should be > 0.\n", canonical_result.id, canonical_result.name, walk_id);
        return 0;
    }

    GibbsSampler *sampler = get_sampler(sampler_id);          //attempt to get the sampler from the map of samplers, if not grab it from the database
//    cout << sampler->to_string() << endl;

    GibbsWalk *walk = new GibbsWalk(walk_id);                 //read the gibbs walk from the database
//    cout << walk->to_string() << endl;

    if (canonical_result.id == 0) {
        //there was some error on this walk
        //update the walk in the database
        ostringstream walk_query;
        walk_query  << "UPDATE gibbs_walk SET "
            << "had_error = true"
            << " WHERE id = " << walk->id;
        mysql_query_check(walk_query.str().c_str());

        cout << "WALK ERRORED!" << endl;
        delete walk;
        return 0;
    }

    //need to see if current samples are in a file
    vector<OUTPUT_FILE_INFO> files;

    int retval = get_output_file_infos(canonical_result, files);
    if (retval) {
        log_messages.printf(MSG_CRITICAL, "[RESULT#%u %s] assimilate handler can't get output filenames for canonical result\n", canonical_result.id, canonical_result.name);  
        return retval;
    }   

    cerr << "files: (" << files.size() << ")" << endl;
    for (uint32_t i = 0; i < files.size(); i++) cerr << "    " << files.at(i).path.c_str() << endl;

    //if there is a samples file, copy it to the data directory
    //file might not exist if there weren't any samples
    if (files.size() == 2 && file_exists(files.at(1).path)) {
        //get the canonical samples file
        //delete all other samples files
        string current_samples_filename = files.at(1).path;

        //copy the samples file to the data directory
        ostringstream saved_samples_filename;
        saved_samples_filename << "/data/dna_at_home/" << sampler->name;

        mkdir(saved_samples_filename.str().c_str(), 0777);  //attempt to make the directory if it doesn't exist
        
        saved_samples_filename << "/walk_" << walk->id << "_steps_" << walk->current_steps;
        cout << "WRITING SAMPLES FILE TO: " << saved_samples_filename.str() << endl;

        ofstream saved_samples_file(saved_samples_filename.str().c_str(), ios::binary);

        ifstream current_samples_file(current_samples_filename.c_str(), ios::binary);

        cout << "current samples file: " << current_samples_filename.c_str() << endl;
//        cout << "contents: " << current_samples_file.rdbuf() << endl;

        saved_samples_file << current_samples_file.rdbuf();

        current_samples_file.close();
        saved_samples_file.close();
    }

//    cout << "CANONICAL RESULT STDERR OUT: " << endl << canonical_result.stderr_out << endl;



    string current_sites;            //read the current sites from the canonical result
    ifstream sites_file(files.at(0).path.c_str());

    string line;
    while (getline(sites_file, line)) {
        current_sites.append(line);
        current_sites.append("\n");
    }   

    //need to remove all windows carriage returns from the file
    current_sites.erase(std::remove(current_sites.begin(), current_sites.end(), '\r'), current_sites.end());

    //the walk hasn't finished, copy the current sites to the download
    //directory for the next workunit in the chain


    //if current sites has a leading newline, get rid of it
    if (current_sites[0] == '\n') {
        current_sites = current_sites.substr(1, current_sites.size() - 2);
    } else {
        current_sites = current_sites.substr(0, current_sites.size() - 1);
    }

//    cout << "current sites: '" << current_sites << "'" << endl;
//    exit(1);

    //update the sites, seed and increment the steps done
    int seed = (int)(drand() * std::numeric_limits<int>::max());
    walk->current_steps += sampler->workunit_steps;
    walk->seed = seed;

    //update the walk in the database
    ostringstream walk_query;
    walk_query  << "UPDATE gibbs_walk SET "
                << "seed = " << seed
                << ", current_steps = " << walk->current_steps
                << ", current_sites = '" << current_sites << "'"
                << " WHERE id = " << walk->id;
    mysql_query_check(walk_query.str().c_str());

    //This walk has finished, don't need to do anything else
    if (walk->current_steps >= sampler->burn_in + sampler->samples) {
        cout << "WALK FINISHED!" << endl;
        delete walk;
        return 0;
    }

    //write sites_filename to disk to be sent with workunit
    ostringstream sites_filename;
    sites_filename << "walk_" << walk->id << "_steps_" << walk->current_steps << "_sites";

    string sites_file_str = sites_filename.str();
    //DONT NEED TO COPY FILE IF WALK HAS FINISHED!!
    char path[1024];
    retval = config.download_path(sites_file_str.c_str(), path);
    cout << "writing sites to: " << path << endl;

    if (retval) return retval;
    FILE* f = fopen(path, "w");
    if (!f) return ERR_FOPEN;
    fprintf(f, "%s", current_sites.c_str());
    fclose(f);

    const char* infiles[2];
    string only_filename = sampler->sequences_filename.substr(sampler->sequences_filename.find_last_of("/\\") + 1);
    infiles[0] = only_filename.c_str();
    infiles[1] = sites_file_str.c_str();

//    cerr << "infile[0]: '" << infiles[0] << "'" << endl;
//    cerr << "infile[1]: '" << infiles[1] << "'" << endl;

    int wu_samples, wu_burn_in;
    //create the next workunit in the walk, if it has not finished
    if (walk->current_steps < sampler->burn_in) {
        //This walk is still doing the burn in

        wu_burn_in = sampler->workunit_steps;
        wu_samples = 0;

        if ((wu_burn_in + walk->current_steps) > sampler->burn_in) {
            //burn in for burn_in = (walk->current_steps + sampler->workunit_steps) - sampler->burn_in
            //sample for samples = (sampler->workunit_steps - burn_in)
            wu_burn_in = (wu_burn_in + walk->current_steps) - sampler->burn_in;   //burn in up to the sampler->burn_in
            wu_samples = (sampler->workunit_steps - wu_burn_in);                 //finish the steps as samples
        }

    } else {
        //this walk has finished the burn in
        wu_burn_in = 0;
        wu_samples = sampler->workunit_steps;
    }

    ostringstream command_line;
    command_line << sampler->command_line_options;
    command_line << " --sample_period " << wu_samples;
    command_line << " --burn_in_period " <<  wu_burn_in;
    command_line << " --current_sites sites.txt ";
    command_line << " --seed " << walk->seed;

//    cout << "the command line is: " << command_line.str() << endl;

    ostringstream additional_xml;
    additional_xml << "<seed>" << seed << "</seed>"
                   << "<sampler_id>" << walk->sampler_id << "</sampler_id>"
                   << "<walk_id>" << walk->id << "</walk_id>"
                   << "<current_steps>" << walk->current_steps << "</current_steps>";

    cerr << "additional_xml: '" << additional_xml.str() << "'" << endl;

    DB_WORKUNIT new_wu;

    ostringstream wu_name;
    wu_name << "gibbs_" << sampler->name << "_" << walk->sampler_id << "_" << "_" << walk->id << "_" << walk->current_steps;

    cout << "creating workunit: '" << wu_name.str() << "'" << endl;

    char name[256];
    sprintf(name, "%s", wu_name.str().c_str());   //should probably make this more descriptive

    // Fill in the job parameters
    new_wu.clear();
    new_wu.appid = wu.appid;
    safe_strcpy(new_wu.name, name);
    new_wu.rsc_fpops_est = sampler->rsc_fpops_est;
    new_wu.rsc_fpops_bound = sampler->rsc_fpops_bound;
    new_wu.rsc_memory_bound = sampler->rsc_memory_bound;
    new_wu.rsc_disk_bound = sampler->rsc_disk_bound; //50MB
    new_wu.delay_bound = sampler->delay_bound;
    new_wu.min_quorum = REPLICATION_FACTOR;
    new_wu.target_nresults = REPLICATION_FACTOR;
    new_wu.max_error_results = REPLICATION_FACTOR*4;
    new_wu.max_total_results = REPLICATION_FACTOR*8;
    new_wu.max_success_results = REPLICATION_FACTOR*4;
    new_wu.batch = sampler_id;


    if (!in_template_initialized) {
        char buf[512];
        sprintf(buf, "templates/%s", sampler->checkpoint_workunit_filename.c_str());
        if (read_file_malloc(config.project_path(buf), in_template)) {
            log_messages.printf(MSG_CRITICAL, "can't read input template %s\n", buf);
            exit(1);
        }
        in_template_initialized = true;
    }


    // Register the job with BOINC
    sprintf(path, "templates/%s", sampler->result_xml_filename.c_str());
    retval = create_work(
        new_wu,
        in_template,                    //input template file
        path,                           //result template file path
        config.project_path(path),
        infiles,                        //input files
        2,                              //number of input files
        config,
        command_line.str().c_str(),
        additional_xml.str().c_str()
    );

    delete walk;

    return 0;
}
