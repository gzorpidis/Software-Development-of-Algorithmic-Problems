#pragma once
#include <iostream>
#include <unistd.h>
#include <string.h>
#include <map>
#include <string>
#include <vector>
#include <assert.h>

//  
//  get_long_opts for Modern C++
//  version 1.0
//  https://github.com/gzorpidis
//
// get_long_options 2023 George Zorpidis <https://gzorpidis.github.io>
// License-Identifier: MIT


/// @brief Parse command line arguments

/*
    @class

    GetLongOption is used to break up (parse) options
    in command lines for easy parsing by shell prodecures
    and more...
    
    
    Usage: 
        (1) instantiate a GetLongOption object,
        passing a vector of strings that are the labels
        of the arguments that are to be used
        e.g. a vector of ["i", "o"] will pick up "-i", "-o"
        in the argument vector (argv)

        Use GetLongOption::options to set how the parsing is done
        SINGLE: Only one argument will be picked up and mapped to 
        each, the very next after the matched label
        MULTIPLE: All arguments will be picked up and mapped to
        each label, until the vector ends or the next matched 
        label is reached.

        Then use parseValues() on the object, passing argc and argv
        And then get_value_from_option or get_vector_from_option to fetch
        the results

        EXAMPLES of "option":

        e.g. `./test -i argv1 -o argv2 argv3`
        labels: [-i, -o]
        if MULTIPLE is set: -o will be mapped/associated to => [argv2, argv3]
        if SINGLE is set, -o will be mapped/associated to => [argv2] and 
        argv3 will be emitted

        for -i, both SINGLE and MULTIPLE will produce the same result,
        mapping -i to [arv1]
    * 
    * 
*/
class Get_long_opts {
    private:

        std::map<std::string, std::string> s_to_s;
        std::map<std::string, std::vector<std::string>> s_to_v;

        bool string_to_string;

        bool isOption(const std::string& option)  {
            return (option[0] == '-' ) ? true : false;
        };

        void addOption (const char* key, const char* value, const int i, const int argc, const char ** argv) {
            if (i + 1 < argc) { 
                if (!isOption(std::string(argv[i+1]))) {
                    if (string_to_string) {
                        s_to_s[std::string(key)] = std::string(value);
                        // std::cout << "inserted: " << key << " -> " << value << std::endl;
                    }
                    else {
                        s_to_v[std::string(value)].push_back(std::string(value));
                    }
                    return;
                }
            } 
            throw std::runtime_error("Option " + std::string(key) + " requires value");
            // fprintf(stderr, "option %s requires value\n", key); exit(EXIT_FAILURE);
        };
    
    public:

        /*
            *   options class should be used in the constructor
            *       to specify SINGLE or MULTIPLE
            *       depending on the desired usage as explained
            *       in the example of the header
        */
        enum class options {
            SINGLE,     // Single item per flag, e.g. -flag yes no , will ignore no
            MULTIPLE    // Multiple item per flag alloed, e.g. -flag yes no will save [yes, no]
        };

        Get_long_opts(std::vector<std::string> labels, options value)  {
            switch (value) {
                case Get_long_opts::options::SINGLE:
                    string_to_string = true;
                    for(const std::string label: labels) {
                        std::string to_insert = "-" + label;
                        s_to_s.insert(std::pair<std::string, std::string>(to_insert, ""));
                    };

                    // for (const auto& pair : s_to_s) {
                    //     std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
                    // }
                    break;
                case Get_long_opts::options::MULTIPLE:
                    string_to_string = false;
                    for(const std::string label: labels) {
                        std::string to_insert = "-" + label;
                        std::vector<std::string> empty_vector;
                        s_to_v[to_insert] = empty_vector;
                    };
                    break;
                default:
                    throw std::runtime_error("Option must be one of SINGLE or MULTIPLE as stated in the Get_long_opts::options enumerator\n");
                    break;
            }
        }
        /*  
        *   @brief Parse values from argc and argv, placing them 
        *   in the calling object for later usage/fetching 
        *   
        *   @param argc     argument counter
        *   @param argv     argument vector
        *   @return         void, parsing the data and placing them in the calling object
        *   
        *   Use parseValues to pass the desired argc and argv
        *   to as parameters to the object, creating the
        *   internal options for later usage
        * 
        */
        void parse_values(int argc, const char** argv) {
    
            for(int i = 0; i < argc; i++) {    
                std::string converted = argv[i];
                if ( isOption(converted) ) {
                    // If string to string is enabled (SINGLE)
                    if (string_to_string) {
                        auto found = s_to_s.find(converted);
                        if (found == s_to_s.end()) {
                            throw std::runtime_error("Invalid option: " + converted);
                        } else {
                            addOption(found->first.c_str(), argv[i+1], i, argc, argv);
                        }
                    } else {
                        // If string to vector is enabled (MULTIPLE)

                        // Make sure it's a legal option
                        auto found = s_to_v.find(converted);
                        if (found == s_to_v.end()) {
                            throw std::runtime_error("Invalid option: " + converted);
                        }

                        // If we run out, throw error
                        int next_token_pos = i + 1;
                        if (next_token_pos >= argc) {
                            throw std::runtime_error("Option " + converted + " requires value");
                        }

                        // If the next option is a option, (non-existente or not) , throw error
                        if (isOption(argv[next_token_pos])) {
                            auto find = s_to_v.find(std::string(argv[next_token_pos]));
                            if (find == s_to_v.end()) {
                                throw std::runtime_error("Invalid option: " + converted);
                            }
                            throw std::runtime_error("Option " + converted + " requires value");
                        }

                        std::string next_token = argv[next_token_pos];

                        while(next_token_pos < argc) {
                            auto found = s_to_v.find(converted);

                            if (found == s_to_v.end()) {
                                throw std::runtime_error("Invalid option: " + converted);
                            } else {
                                s_to_v[std::string(converted)].push_back(argv[next_token_pos]);
                            }

                            next_token_pos++;

                            if (next_token_pos >= argc) break;

                            next_token = argv[next_token_pos];

                            if (isOption(next_token)) {break;}
                        }
                    }
                }
            }
        }
        




        /*
        *   @brief Return a single value mapped to the option/label provided as a parameter
        *   
        * 
        *   @param option The option you want to find the value of
        *   @return A single <const char*> of the associated value of the option,as it has been parsed, use delete after usage,
        *   or NULL if option is not part of options as set in the constructor, or if the option has no associated values
        *   
        * 
        *   Given an option value, search and return for its associated value(s).
        *   If the label/option is associated to many values, the first one is returned, as it has been parsed in linear order
        *   If the label/option is associated to no values, NULL is returned
        *   If the label/option is not part of the labels vector given in the constructor, NULL is returned
        */
        char* get_value_from_option(const char* option)  {
          std::string search;

          //  option must not depend on starting "-"
          //      if it does not exist, add the hyphen "-"
          //      or do not add anything if it already exists
          if (!isOption(option)) {
              search = "-" + std::string(option);
          } else {
              search = option;
          }

          if (string_to_string) {
              auto found = s_to_s.find(search);
              if (found != s_to_s.end()) {
                  std::string output = found->second;
                  if (output.empty()) {
                      return NULL;
                  }
                  char* to_return = new char[output.size() + 1];
                  strncpy((to_return), output.c_str(), output.size()+1);
                  return to_return;
              } else {
                  return NULL;
              }
          } else {
              auto found = s_to_v.find(search);
              if (found != s_to_v.end()) {
                  if (found->second.empty()) {
                      return NULL;
                  } else {
                      std::string output = found->second[0];
                      char* to_return = new char[output.size() + 1];
                      strncpy((to_return), output.c_str(), output.size()+1);
                      return to_return;
                  }
              } else {
                  return NULL;
              }
          }

          return NULL;
        }       



        /*
        *   *@brief Return a vector of values mapped to the option/label provided as a parameter
        *   
        *   Given an option value, search and return for its associated value(s).
        *   If the label/option is associated to many values, all of them are returned in linear order as they have been parsed
        *   If the label/option is associated to no values, empty vector is returned
        *   If the label/option is not part of the labels vector given in the constructor, empty vector is returned
        *   @param option The option you want to find the value of
        *
        *   @return A vector of <std::string> of the associated values of the option,as it has been parsed, use delete after usage,
        *   or an empty vector if option is not part of options as set in the constructor, or if the option has no associated values
        *   
        */
        std::vector<std::string> get_vector_from_option(const char* option)  {

            std::vector<std::string> to_return;
            std::string search;

            if (!isOption(option)) {
                search = "-" + std::string(option);
            } else {
                search = option;
            }

            if (string_to_string) {
                const char* single = get_value_from_option(search.c_str());

                if (!single) return to_return;
                else {
                    to_return.push_back(single);
                    delete[] single;
                    return to_return;
                }
            } else {
                auto find = s_to_v.find(search);

                if (find == s_to_v.end()) {
                    return to_return;
                } else {
                    for(const std::string& value : find->second ) {
                        to_return.push_back(value);
                    }
                    return to_return;
                }
            }

            assert(to_return.empty());
            return to_return;
        }

        /*
        *   @brief Print the state of the associations/mappings of the labels and options
        *   
        *   @return Prints the mappings in the console, useful for debugging
        * 
        
        */
        void print() {
            if (!string_to_string) {
                for (const auto& pair : s_to_v) {
                    std::cout << "Key: " << pair.first << ", Values: ";
                    for(const auto& item: pair.second) {
                        std::cout << item << " ";
                    }
                    std::cout << std::endl;
                }
            } else {
                for(const auto& pair : s_to_s) {
                    std::cout << "Key: " << pair.first << ", Value -> " << pair.second << std::endl;
                }
                
            }
        }
};