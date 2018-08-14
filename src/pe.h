#include <iostream>
#include <string>
#ifndef PRESENTATION_ERROR_H_
#define PRESENTATION_ERROR_H_
class pe:public std::logic_error
{
public:
	explicit pe(const std::string& what_arg):std::logic_error(what_arg){};
 	explicit pe(const char* what_arg):std::logic_error(what_arg){};
};
#endif