#ifndef MEM_EKF_EXAMPLE_CSV_READER_HPP_
#define MEM_EKF_EXAMPLE_CSV_READER_HPP_

#include <fstream>
#include <sstream>
#include <vector>
#include <string>

class CsvReader {
  public:
    CsvReader(std::string file_name, std::string delimeter = std::string(","))
      : file_name_(file_name)
      , delimeter_(delimeter) {}

    std::vector<std::vector<std::string>> GetData(void);
    
  private:
    std::string file_name_;
    std::string delimeter_;
};


#endif  //  MEM_EKF_EXAMPLE_CSV_READER_HPP_
