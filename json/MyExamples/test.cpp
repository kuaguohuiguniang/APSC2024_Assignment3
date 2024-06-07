#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "json.hpp"
int main()
{
  // The json object type
  using json=nlohmann::json;
  // Read from a file
  std::ifstream ifile("data.json");
  json j1;
  ifile >>j1;
  // print on screen (use setw(4) for pretty print)
  std::cout<<std::setw(4)<<j1<<std::endl;
  // Extract data. you need to use get<type> to select the type
  // Strings
  auto s = j1["name"].get<std::string>();
  std::cout<<s<<std::endl;
  // vectors
  auto v = j1["list"].get<std::vector<double>>();
  auto vs =j1["stringlist"].get<std::vector<std::string>>();
  for (auto x : v)std::cout<<x<<" ";
  std::cout<<std::endl;
  for (auto x : vs)std::cout<<x<<" ";
 std::cout<<std::endl;
  //boolean
//  auto happy=j1["happy"].get<bool>();
  // Alternative with default
  auto happy=j1.value("happy",false);
  std::cout<<std::boolalpha<<happy<<std::endl;
  auto sad=j1.value("sad",false);

  std::cout<<std::boolalpha<<sad<<std::endl;
  
    //more complex structures

  auto money = j1["object"]["value"].get<double>();
  std::cout<<money<<std::endl;
// write on file (compact form)
  std::ofstream ofile("output.json");
  ofile<<j1;
  // pretty output
  ofile.close();
  ofile.open("output_pretty.json");
  ofile<<std::setw(4)<<j1;
  
  
}
