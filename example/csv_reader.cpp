#include "csv_reader.hpp"



std::vector<std::vector<std::string>> CsvReader::GetData(void) {
    std::ifstream file(file_name_);
    std::vector<std::vector<std::string>> data_list;
    std::string line = std::string("");

    // Iterate through each line and split the content using delimeter
    while (getline(file, line)) {
      std::vector<std::string> vec;
      std::stringstream ss(line);

      while (ss.good()) {
          std::string substr;
          getline(ss, substr, delimeter_[0]);
          vec.push_back(substr);
      }

      data_list.push_back(vec);
    }

    // Close the File
    file.close();

    return data_list;
}
