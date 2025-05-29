/**
@file Exceptions.h
@brief Declaration of the various exception classes.

THAMES tries to implement exception handling consistently throughout the code,
although more could be done to check for out-of-bounds errors on arrays and
other container classes.

Right now, many of these classes are quite similar and could probably be
structured better with a generic base class, but the thought is that each class
will take on more distinct details later.

@todo Make a virtual base class for common things like descriptions, etc.

*/

#ifndef SRC_THAMESLIB_EXCEPTIONS_H_
#define SRC_THAMESLIB_EXCEPTIONS_H_

#include <iostream>
#include <string>

/**
@class Declare the EOBException class

The `EOBException` class handles exceptions caused by attempting to access
an out of bounds element of an array.
*/
class EOBException {

private:
  std::string arrayname_;    /**< Name of the array accessed */
  std::string classname_;    /**< Name of the class that accessed the array */
  std::string functionname_; /**< Name of the method that accessed the array */
  int sizelimit_;            /**< Number of elements contained in the array */
  int indx_; /**< Out-of-bounds element number that was queried */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  EOBException() {
    classname_ = "";
    functionname_ = "";
    arrayname_ = "";
    sizelimit_ = 0;
    indx_ = 0;
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param arname is the name of the array that was accessed
  @param sl is the total number of array elements in that array
  @param id is the element number (out of bounds) that was queried erroneously
  */
  EOBException(const std::string &cname, const std::string &fileName,
               const std::string &arname, const int sl, const unsigned int id) {
    classname_ = cname;
    functionname_ = fileName;
    arrayname_ = arname;
    sizelimit_ = sl;
    indx_ = id;
  }

  /**
  @brief Get the class name responsible for throwing the EOB exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the EOB exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the array name that was queried when the EOB exception was thrown.

  @return the array name
  */
  std::string &getArrayname() const { return (std::string &)arrayname_; }

  /**
  @brief Get the number of elements of the array queried when the EOB exception
  was thrown.

  @return the array name
  */
  int getSizelimit() const { return sizelimit_; }

  /**
  @brief Get the queried (out-of-bounds) index number.

  @return the queried index number
  */
  unsigned int getIndx() const { return indx_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "EOB Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << std::endl << "EOB Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    if (indx_ == 0) {
      std::cout << "        Array: " << arrayname_ << std::endl;
      std::cerr << "        Array: " << arrayname_ << std::endl;
    } else {
      std::cout << "        Array: " << arrayname_ << " contains "
                << sizelimit_;
      std::cout << " elements, but tried to access element " << indx_
                << std::endl;
      std::cerr << "        Array: " << arrayname_ << " contains "
                << sizelimit_;
      std::cerr << " elements, but tried to access element " << indx_
                << std::endl;
    }
    return;
  }

}; // End of the EOBException class

/**
@class Declare the FileException class

The `FileException` class handles exceptions caused by attempting to open,
close, write to, or read from a file.

*/
class FileException {

private:
  std::string filename_;     /**< Name of the offending file */
  std::string extype_;       /**< Name of the exception description */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  FileException() {
    classname_ = "";
    functionname_ = "";
    filename_ = "";
    extype_ = "";
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param filename is the name of the offending file
  @param extype is the description of the exception type
  */
  FileException(const std::string &cname, const std::string &fileName,
                const std::string &filename, const std::string &extype) {
    classname_ = cname;
    functionname_ = fileName;
    filename_ = filename;
    extype_ = extype;
  }

  /**
  @brief Get the class name responsible for throwing the file exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the file exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the file name that was queried when the file exception was thrown.

  @return the file name
  */
  std::string &getFilename() const { return (std::string &)filename_; }

  /**
  @brief Get the file exception type description.

  @return the exception type description
  */
  std::string &getExtype() const { return (std::string &)extype_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "File Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        File: " << filename_ << ", Problem:" << extype_
              << std::endl;
    std::cout << std::endl << "File Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cerr << "        File: " << filename_ << ", Problem: " << extype_
              << std::endl;
    return;
  }

}; // End of the FileException class

/**
@class Declare the FloatException class

The `FloatException` class handles exceptions caused by floating point
operations, especially divide-by-zero exceptions.

*/
class FloatException {

private:
  std::string description_;  /**< Description of the floating point exception */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  FloatException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param strd is the description of the exception
  */
  FloatException(const std::string &cname, const std::string &fileName,
                 const std::string &strd) {
    classname_ = cname;
    functionname_ = fileName;
    description_ = strd;
  }

  /**
  @brief Get the class name responsible for throwing the floating point
  exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the floating point
  exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the description of the floating point exception.

  @return the file name
  */
  std::string &getDescription() const { return (std::string &)description_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "Floating Point Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        Description: " << description_ << std::endl;
    std::cerr << std::endl << "Floating Point Exception Thrown:" << std::endl;
    std::cerr << "    Details: " << std::endl;
    std::cerr << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cerr << "        Description: " << description_ << std::endl;
    return;
  }

}; // End of the FloatException class

/**
@class Declare the HandleException class

The `HandleException` class handles exceptions caused by errors in dealing
with data handles.

*/
class HandleException {

private:
  std::string description_;  /**< Description of the handle exception */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */
  std::string handle_; /**< Description of the handle causing the exception */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  HandleException() {
    classname_ = "";
    functionname_ = "";
    handle_ = "";
    description_ = "";
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param handle is the handle that caused the exception
  @param strd is the description of the exception
  */
  HandleException(const std::string &cname, const std::string &fileName,
                  const std::string &handle, const std::string &strd) {
    classname_ = cname;
    functionname_ = fileName;
    handle_ = handle;
    description_ = strd;
  }

  /**
  @brief Get the class name responsible for throwing the handle exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the handle exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the handle causing the exception.

  @return the handle
  */
  std::string &getHandle() const { return (std::string &)handle_; }

  /**
  @brief Get the description of the handle exception.

  @return the file name
  */
  std::string &getDescription() const { return (std::string &)description_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "Handle Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        Description: " << description_ << std::endl;
    std::cout << "             Handle: " << handle_ << std::endl;
    std::cout << std::endl << "Floating Point Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        Description: " << description_ << std::endl;
    std::cout << "             Handle: " << handle_ << std::endl;
    return;
  }

}; // End of the HandleException class

/**
@class Declare the GEMException class

The `GEMException` class handles exceptions caused by errors originating
in the GEM3K library.

*/
class GEMException {

private:
  std::string description_;  /**< Description of the GEM exception */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  GEMException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param strd is the description of the exception
  */
  GEMException(const std::string &cname, const std::string &fileName,
               const std::string &strd) {
    classname_ = cname;
    functionname_ = fileName;
    description_ = strd;
  }

  /**
  @brief Get the class name responsible for throwing the handle exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the handle exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the description of the handle exception.

  @return the file name
  */
  std::string &getDescription() const { return (std::string &)description_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "GEM Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        " << description_ << std::endl;
    std::cout << std::endl << "GEM Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        " << description_ << std::endl;
    return;
  }

}; // End of GEMException class

/**
@class Declare the MicrostructureException class

The `MicrostructureException` class handles exceptions related to
microstructure modifications.

*/
class MicrostructureException {
private:
  std::string description_;  /**< Description of the GEM exception */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */
  bool excp_; /**< true <-> for exception / false <-> for normal exit */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.
*/
  MicrostructureException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
    excp_ = true;
  }

  /**
    @brief Overloaded constructor that is typically invoked by THAMES.

    @param cname is the class name where the exception was thrown
    @param fileName is the method name where the exception was thrown
    @param strd is the description of the exception
  */
  MicrostructureException(const std::string &cname, const std::string &fileName,
                          const std::string &strd) {
    classname_ = cname;
    functionname_ = fileName;
    description_ = strd;
  }

  /**
    @brief Overloaded constructor that is typically invoked by THAMES.

    @param cname is the class name where the exception was thrown
    @param fileName is the method name where the exception was thrown
    @param strd is the description of the exception
    @param excp is true <-> for exception / false <-> for normal exit
  */
  MicrostructureException(const std::string &cname, const std::string &fileName,
                          const std::string &strd, bool excp) {
    classname_ = cname;
    functionname_ = fileName;
    description_ = strd;
    excp_ = excp;
  }

  /**
  @brief Get the class name responsible for throwing the exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the description of the exception.

  @return the file name
  */
  std::string &getDescription() const { return (std::string &)description_; }

  bool getExcp() { return excp_; }

  /**
@brief Provide formatted output of the exception details.

*/
  void printException() {
    // bool excp1_ = true;
    if (excp_) {
      std::cout << std::endl << "Microstructure Exception Thrown:" << std::endl;
      std::cout << "    Details: " << std::endl;
      std::cout << "        Offending Function " << classname_
                << "::" << functionname_ << std::endl;
      std::cout << "        Problem: " << description_ << std::endl;
      std::cout << std::endl << "Microstructure Exception Thrown:" << std::endl;
      std::cout << "    Details: " << std::endl;
      std::cout << "        Offending Function " << classname_
                << "::" << functionname_ << std::endl;
      std::cout << "        Problem: " << description_ << std::endl;
    } else {
      std::cout << std::endl << "Microstructure Exception Thrown:" << std::endl;
      std::cout << "    Details: " << std::endl;
      std::cout << "        From Function " << classname_
                << "::" << functionname_ << std::endl;
      std::cout << "        reason: " << description_ << std::endl;
      std::cout << std::endl << "Microstructure Exception Thrown:" << std::endl;
      std::cout << "    Details: " << std::endl;
      std::cout << "        From Function " << classname_
                << "::" << functionname_ << std::endl;
      std::cout << "        reason: " << description_ << std::endl;
    }

    return;
  }
}; // End of MicrostructureException class

/**
@class Declare the DataException class

The `DataException` class handles exceptions related to miscellaneous
data errors.

*/
class DataException {

private:
  std::string description_;  /**< Description of the GEM exception */
  std::string classname_;    /**< Name of the class that threw the exception */
  std::string functionname_; /**< Number of function that threw the exception */

public:
  /**
  @brief Default constructor initializes class members to default (blank)
  values.

  */
  DataException() {
    classname_ = "";
    functionname_ = "";
    description_ = "";
  }

  /**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param strd is the description of the exception
  */
  DataException(const std::string &cname, const std::string &functionName,
                const std::string &strd) {
    classname_ = cname;
    functionname_ = functionName;
    description_ = strd;
  }

  /**
  @brief Get the class name responsible for throwing the handle exception.

  @return the class name
  */
  std::string &getClassname() const { return (std::string &)classname_; }

  /**
  @brief Get the function name responsible for throwing the handle exception.

  @return the function name
  */
  std::string &getFunctionname() const { return (std::string &)functionname_; }

  /**
  @brief Get the description of the handle exception.

  @return the file name
  */
  std::string &getDescription() const { return (std::string &)description_; }

  /**
  @brief Provide formatted output of the exception details.

  */
  void printException() {
    std::cout << std::endl << "Data Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        Problem:" << description_ << std::endl;
    std::cout << std::endl << "Data Exception Thrown:" << std::endl;
    std::cout << "    Details: " << std::endl;
    std::cout << "        Offending Function " << classname_
              << "::" << functionname_ << std::endl;
    std::cout << "        Problem: " << description_ << std::endl;
    return;
  }

}; // End of DataException class

#endif // SRC_THAMESLIB_EXCEPTIONS_H_
