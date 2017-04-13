#ifndef NOTYETIMPLEMENTEDERROR_H
#define NOTYETIMPLEMENTEDERROR_H

#include <stdexcept>
#include <string>

namespace ctbn {

    using std::string;
    using std::logic_error;

    class not_yet_implemented_error : logic_error {
    public:
        not_yet_implemented_error(
         string const & the_class_name,
         string const & the_method_name) :
        logic_error(string("The function ") +
         the_class_name + (the_class_name == "" ? "" : "::") +
         the_method_name + " has not yet been implemented.") {};
    };
}

#endif
