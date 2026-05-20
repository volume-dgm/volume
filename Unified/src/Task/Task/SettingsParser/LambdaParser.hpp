#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <functional>
#include <memory>
#include <cmath>

struct Thunk;

class GenericDataType;

GenericDataType HelperFunc(std::shared_ptr<Thunk> v);
GenericDataType HelperFunc2(std::shared_ptr<Thunk> v, bool log_debug = false, int layer = 0);
std::vector<std::string> splitString(const std::string& input, const std::string& delimiters, const std::string& delimitersToOutput) {
    std::vector<std::string> tokens;
    std::size_t startPos = 0;
    std::size_t endPos;

    const std::string delimitersAll = delimiters + delimitersToOutput;

    while ((endPos = input.find_first_of(delimitersAll, startPos)) != std::string::npos) {
        // Extract substring between startPos and the found delimiter
        if (endPos > startPos) {
            tokens.push_back(input.substr(startPos, endPos - startPos));
        }
        if (endPos < input.length()) {
          if (delimitersToOutput.find(input[endPos]) != std::string::npos) {
            // should parse delimiter as own thing
            tokens.push_back(std::string() + input[endPos]); // conversion of char to string, probably not the best one
          }
        }
        // Advance start position past the delimiter
        startPos = endPos + 1;
    }

    // Add the last part of the string if it doesn't end with a delimiter
    if (startPos < input.length()) {
        tokens.push_back(input.substr(startPos));
    }

    return tokens;
}

class GenericDataType {
  public:
    enum DataType {
      NUMBER,
      STRING,
      BOOL,
      NIL,
      FUNCTION,
      HASHMAP,
      VARIABLE, // intermediate type for incomplete lambda expressions
      FUNCTIONAPPL // intermediate type for incomplete lambda expressions
    };

    struct GenericDataTypeHash {
      std::size_t operator()(const GenericDataType& k) const noexcept {
        if (k.myType == HASHMAP) { return 0; } // TOCHANGE?
        if (k.myType == FUNCTION) { return 0; } // TOCHANGE?
        if (k.myType == NUMBER) { return std::hash<double>{}(k.num_val); }
        if (k.myType == STRING) { return std::hash<std::string>{}(k.str_val); }
        if (k.myType == BOOL) { return std::hash<bool>{}(k.bool_val); }
        if (k.myType == NIL) { return 1; }
        if (k.myType == VARIABLE) { return std::hash<std::string>{}(k.var_val); }
        if (k.myType == FUNCTIONAPPL) {
            GenericDataTypeHash tmp = GenericDataTypeHash();
            // return 0;
            return tmp(HelperFunc((*k.funcAppl_val).first)) ^ (tmp(HelperFunc((*k.funcAppl_val).second)) << 1);
        }
        //throw std::runtime_error("unknown type in hash");
        return 1;
      }
    };

    using HashMap = std::unordered_map<GenericDataType, GenericDataType, GenericDataTypeHash>;
    using FuncType = std::function<GenericDataType(std::shared_ptr<Thunk>)>;
    using FuncTypeAppl = std::pair<std::shared_ptr<Thunk>, std::shared_ptr<Thunk>>;

    DataType myType = NIL;
    double num_val = 0;
    std::string str_val = "";
    bool bool_val = false;
    std::shared_ptr<FuncType> func_val;
    std::shared_ptr<HashMap> hashmap_val = {};
    std::shared_ptr<FuncTypeAppl> funcAppl_val;
    std::string var_val = "";
    std::string func_comment = "";

    bool operator==(const GenericDataType& other) const {
      if (myType != other.myType) { return false; }
      if (myType == HASHMAP) { return hashmap_val == other.hashmap_val; }
      if (myType == FUNCTION) { return func_val == other.func_val; }
      if (myType == NUMBER) { return num_val == other.num_val; }
      if (myType == STRING) { return str_val == other.str_val; }
      if (myType == BOOL) { return bool_val == other.bool_val; }
      if (myType == NIL) { return true; }
      if (myType == VARIABLE) { return var_val == other.var_val; }
      if (myType == FUNCTIONAPPL) { return funcAppl_val == other.funcAppl_val; }
      throw std::runtime_error("unknown type in eqaulity");
    }

    GenericDataType() {}

    GenericDataType(std::string s) {
      myType = STRING;
      str_val = s;
    }

    GenericDataType(double d) {
      myType = NUMBER;
      num_val = d;
    }

    /*GenericDataType(bool b) {
      myType = BOOL;
      bool_val = b;
    }*/


    std::string displayStr() {
      if (myType == HASHMAP) {
          std::string res = "";
          for (auto it = (*hashmap_val).begin(); it != (*hashmap_val).end(); ++it) {
            GenericDataType f = it -> first;
            GenericDataType s = it -> second;
            res += "KEY: " + f.displayStr() + " VALUE: " + s.displayStr() + ", ";
          }
          return "HASHMAP: [" + res + "]";
      }
      if (myType == FUNCTION) { return "FUNCTION: " + func_comment; }
      if (myType == NUMBER) { return "NUMBER: " + std::to_string(num_val); }
      if (myType == STRING) { return "STRING: " + str_val; }
      if (myType == BOOL) { return "BOOL: " + std::to_string(bool_val); }
      if (myType == NIL) { return "NIL"; }
      if (myType == VARIABLE) { return "VARIABLE: " + var_val; }
      if (myType == FUNCTIONAPPL) { return "FUNCTION APPLICATION: (" + HelperFunc((*funcAppl_val).first).displayStr() + ") APPLIED TO (" + HelperFunc((*funcAppl_val).second).displayStr() + ")"; }
      throw std::runtime_error("unknown type in display string");
    }

    GenericDataType eval(bool log_debug = false, int inner_layer = 0) {
      if (myType == HASHMAP ||
          myType == FUNCTION ||
          myType == NUMBER ||
          myType == STRING ||
          myType == BOOL ||
          myType == NIL) {
            return *this;
          }
      if (myType == VARIABLE) { throw std::runtime_error("tried to evaluate unbound variable: " + var_val); }
      if (myType == FUNCTIONAPPL) {
        GenericDataType fun = HelperFunc2((*funcAppl_val).first, log_debug, inner_layer + 1);
        if (log_debug) {
          std::cout << "[EVAL LOG] " + std::string(inner_layer, '|') + "the function is evaluated as: " + fun.displayStr() << std::endl;
        }
        if (fun.myType != FUNCTION) {
          throw std::runtime_error("tried to call non-function object: " + fun.displayStr());
        }
        auto arg = (*funcAppl_val).second; // eval must be called inside the function if needed
        // std::cout << "calling function with argument: " + arg.displayStr() << std::endl;
        // (*fun.func_val)(arg);
        // std::cout << "call successful" << std::endl;
        // std::cout << "res is: " + (*fun.func_val)(arg).displayStr() << std::endl;

        return (*fun.func_val)(arg).eval(log_debug, inner_layer + 1);
        //return (*fun.func_val)(arg).eval(log_debug, inner_layer + 1);
      }
      throw std::runtime_error("unknown type in eval");
    }
};



struct Thunk {
  GenericDataType expr;
  bool evaluated = false;
  GenericDataType value;

  Thunk(GenericDataType expr_in) {
    expr = expr_in;
  }

  GenericDataType force(bool log_debug = false, int layer = 0) {
    if (!evaluated) {
      value = expr.eval(log_debug, layer + 1);
      evaluated = true;
    }
    return value;
  }
};

GenericDataType HelperFunc(std::shared_ptr<Thunk> v) {
  return v->expr;
}
GenericDataType HelperFunc2(std::shared_ptr<Thunk> v, bool log_debug, int layer) {
  return v->force(log_debug, layer);
}
class FunctionEntry {
  public:
    using FuncType = std::function<GenericDataType(std::shared_ptr<Thunk>)>;
    using FuncTypeAppl = std::pair<std::shared_ptr<Thunk>, std::shared_ptr<Thunk>>;
    std::string description = "";
    std::string name = "";
    std::shared_ptr<FuncType> func_val;
    FunctionEntry(std::string func_name, std::shared_ptr<FuncType> func, std::string desc) {
      description = desc;
      name = func_name;
      func_val = func;
    }
    FunctionEntry(std::string func_name, FuncType func, std::string desc) {
      description = desc;
      name = func_name;
      func_val = std::make_shared<FuncType>(func);
    }


    static std::shared_ptr<FuncType> getCurriedFunction(std::shared_ptr<std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)>> original) {
      auto fun = ([original](std::shared_ptr<Thunk> x) -> GenericDataType {
        auto x_val = x;
        auto fun = [original, x_val](std::shared_ptr<Thunk> y) -> GenericDataType {
          return (*original)(x_val, y);
        };
        GenericDataType tmp = GenericDataType();
        tmp.myType = GenericDataType::FUNCTION;
        tmp.func_val = std::make_shared<FuncType>(fun);
        return tmp;
      });
      return std::make_shared<FuncType>(fun);
    }

    static std::shared_ptr<FuncType> getCurriedFunction(std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)> original) {
      auto fun = std::make_shared<std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)>>(original);
      return getCurriedFunction(fun);
    }

    static std::shared_ptr<FuncType> getCurriedFunctionThree(std::shared_ptr<std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)>> original) {
      auto fun = ([original](std::shared_ptr<Thunk> x) -> GenericDataType {
        auto x_val = x;
        auto fun = getCurriedFunction([original, x_val](std::shared_ptr<Thunk> y, std::shared_ptr<Thunk> z) -> GenericDataType {
          return (*original)(x_val, y, z);
        });
        GenericDataType tmp = GenericDataType();
        tmp.myType = GenericDataType::FUNCTION;
        tmp.func_val = fun;
        return tmp;
      });
      return std::make_shared<FuncType>(fun);
    }

    static std::shared_ptr<FuncType> getCurriedFunctionThree(std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)> original) {
      auto fun = std::make_shared<std::function<GenericDataType(std::shared_ptr<Thunk>, std::shared_ptr<Thunk>, std::shared_ptr<Thunk>)>>(original);
      return getCurriedFunctionThree(fun);
    }

    static std::vector<FunctionEntry> getDefaultFunctions() {
      std::vector<FunctionEntry> res;
      // ARITHMETICS
      res.push_back(FunctionEntry(
        "+",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            return GenericDataType(xval.num_val + yval.num_val);
          }
          throw std::runtime_error("+ must be called with numbers");
        }),
        "adds together two numbers"
        ));
      res.push_back(FunctionEntry(
        "-",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            return GenericDataType(xval.num_val - yval.num_val);
          }
          throw std::runtime_error("- must be called with numbers");
        }),
        "subtracts two numbers"
        ));
      res.push_back(FunctionEntry(
        "*",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            return GenericDataType(xval.num_val * yval.num_val);
          }
          throw std::runtime_error("* must be called with numbers");
        }),
        "multiplies two numbers"
        ));
      res.push_back(FunctionEntry(
        "/",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            return GenericDataType(xval.num_val / yval.num_val);
          }
          throw std::runtime_error("/ must be called with numbers");
        }),
        "divides two numbers"
        ));
      res.push_back(FunctionEntry(
        "sqrt",
        [](std::shared_ptr<Thunk> x) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::NUMBER) {
            return GenericDataType(std::sqrt(xval.num_val));
          }
          throw std::runtime_error("sqrt must be called with number");
        },
        "takes square root of number"
        ));
      // COMPARE
      res.push_back(FunctionEntry(
        ">",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.num_val > yval.num_val;
            return res;
          }
          throw std::runtime_error("> must be called with numbers");
        }),
        "compares two numbers"
        ));
      res.push_back(FunctionEntry(
        "<",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.num_val < yval.num_val;
            return res;
          }
          throw std::runtime_error("< must be called with numbers");
        }),
        "compares two numbers"
        ));
      res.push_back(FunctionEntry(
        ">=",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.num_val >= yval.num_val;
            return res;
          }
          throw std::runtime_error(">= must be called with numbers");
        }),
        "compares two numbers"
        ));
      res.push_back(FunctionEntry(
        "<=",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::NUMBER && yval.myType == GenericDataType::NUMBER) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.num_val <= yval.num_val;
            return res;
          }
          throw std::runtime_error("<= must be called with numbers");
        }),
        "compares two numbers"
        ));
      res.push_back(FunctionEntry(
        "==",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          GenericDataType res = GenericDataType();
          res.myType = GenericDataType::BOOL;
          res.bool_val = (xval == yval);
          return res;
        }),
        "compares two values"
        ));
      // BOOL
      res.push_back(FunctionEntry(
        "and",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::BOOL && yval.myType == GenericDataType::BOOL) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.bool_val && yval.bool_val;
            return res;
          }
          throw std::runtime_error("and must be called with bools");
        }),
        "logical and"
        ));
      res.push_back(FunctionEntry(
        "or",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::BOOL && yval.myType == GenericDataType::BOOL) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.bool_val || yval.bool_val;
            return res;
          }
          throw std::runtime_error("or must be called with bools");
        }),
        "logical or"
        ));
      res.push_back(FunctionEntry(
        "xor",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          GenericDataType yval = y -> force();
          if (xval.myType == GenericDataType::BOOL && yval.myType == GenericDataType::BOOL) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = xval.bool_val ^ yval.bool_val;
            return res;
          }
          throw std::runtime_error("xor must be called with bools");
        }),
        "logical xor"
        ));
      res.push_back(FunctionEntry(
        "not",
        [](std::shared_ptr<Thunk> x) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::BOOL) {
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::BOOL;
            res.bool_val = !xval.bool_val;
            return res;
          }
          throw std::runtime_error("not must be called with bools");
        },
        "logical not"
        ));
      // HASHMAP
      res.push_back(FunctionEntry(
        "insert",
        getCurriedFunctionThree([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y, std::shared_ptr<Thunk> z) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::HASHMAP) {
            GenericDataType yval = y -> force();
            GenericDataType zval = z -> force();
            auto newMap = (*xval.hashmap_val);
            newMap[yval] = zval;
            GenericDataType res = GenericDataType();
            res.myType = GenericDataType::HASHMAP;
            res.hashmap_val =  std::make_shared<std::unordered_map<GenericDataType, GenericDataType, GenericDataType::GenericDataTypeHash>>(newMap);
            return res;
          }
          throw std::runtime_error("first argument of insert must be a hashmap, instead was " + xval.displayStr());
        }),
        "insert hashmap key val"
        ));
      res.push_back(FunctionEntry(
        "get",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::HASHMAP) {
            GenericDataType yval = y -> force();
            return (*xval.hashmap_val)[yval];
          }
          throw std::runtime_error("first argument of get must be a hashmap, instead was " + xval.displayStr());
        }),
        "get hashmap key"
        ));
      // META
      res.push_back(FunctionEntry(
        "type",
        [](std::shared_ptr<Thunk> x) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::HASHMAP) { return GenericDataType("HASHMAP"); }
          if (xval.myType == GenericDataType::FUNCTION) { return GenericDataType("FUNCTION"); }
          if (xval.myType == GenericDataType::NUMBER) { return GenericDataType("NUMBER"); }
          if (xval.myType == GenericDataType::STRING) { return GenericDataType("STRING"); }
          if (xval.myType == GenericDataType::BOOL) { return GenericDataType("BOOL"); }
          if (xval.myType == GenericDataType::NIL) { return GenericDataType("NIL"); }
          if (xval.myType == GenericDataType::VARIABLE) { return GenericDataType("VARIABLE"); }
          if (xval.myType == GenericDataType::FUNCTIONAPPL) { return GenericDataType("FUNCTION APPLICATION"); }
          throw std::runtime_error("unknown type in type");
        },
        "outputs the type of the argument as string"
        ));
      res.push_back(FunctionEntry(
        "ignore",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          return y -> force();
        }),
        "ignores the first argument, useful as a comment, same as lambda x y.y"
        ));
      res.push_back(FunctionEntry(
        "ignoreDebug",
        getCurriedFunction([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y) -> GenericDataType {
          x -> force();
          return y -> force();
        }),
        "ignores the first argument, but evaluates it, useful as a way to debug"
        ));
      res.push_back(FunctionEntry(
        "debug",
        [](std::shared_ptr<Thunk>x) -> GenericDataType {
          GenericDataType xval = x -> force();
          std::cout << "[DEBUG] " << xval.displayStr() << std::endl;
          return xval;
        },
        "outputs the input as string to console, outputting the input itself"
        ));
      // CONTROL
      res.push_back(FunctionEntry(
        "if",
        getCurriedFunctionThree([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y, std::shared_ptr<Thunk> z) -> GenericDataType {
          GenericDataType xval = x -> force();
          if (xval.myType == GenericDataType::BOOL) {
            if (xval.bool_val) {
              return y -> force();
            } else {
              return z -> force();
            }
          }
          throw std::runtime_error("first argument of if must be bool, instead was " + xval.displayStr());
        }),
        "if cond executeIfTrue executeIfFalse"
        ));
      res.push_back(getHelpFunction(res));
      res.push_back(FunctionEntry(
        "while",
        getCurriedFunctionThree([](std::shared_ptr<Thunk> x, std::shared_ptr<Thunk> y, std::shared_ptr<Thunk> z) -> GenericDataType {
          GenericDataType cond = GenericDataType();
          cond.myType = GenericDataType::FUNCTIONAPPL;
          std::shared_ptr<Thunk> z_copy = std::make_shared<Thunk>(z -> force());
          FuncTypeAppl res = FuncTypeAppl(x, z_copy);
          cond.funcAppl_val = std::make_shared<FuncTypeAppl>(res); // apply the function to the argument
          while (true) {
            // std::cout << "while looped checking" << std::endl;
            GenericDataType tmp_cond = cond.eval();
            if (tmp_cond.myType != GenericDataType::BOOL) {
              throw std::runtime_error("condigiton argument of if must be bool, instead was " + tmp_cond.displayStr());
            }
            if (!tmp_cond.bool_val) {
              // std::cout << "while looped ending" << std::endl;
              break;
            }
            // std::cout << "while looped looping" << std::endl;
            GenericDataType tmp = GenericDataType();
            tmp.myType = GenericDataType::FUNCTIONAPPL;
            FuncTypeAppl res2 = FuncTypeAppl(y, z_copy);
            tmp.funcAppl_val = std::make_shared<FuncTypeAppl>(res2); // apply the function to the argument
            z_copy = std::make_shared<Thunk>(tmp.eval());
            FuncTypeAppl res3 = FuncTypeAppl(x, z_copy);
            cond.funcAppl_val = std::make_shared<FuncTypeAppl>(res3); // apply the function to the argument
          }
          return z_copy -> force();
        }),
        "while cond body state\nevaluates cond(state), returns state if false, else changes the state to be (body state) and runs again"
        ));
      res.push_back(getHelpFunction(res));
      return res;
    }

    static FunctionEntry getHelpFunction(std::vector<FunctionEntry>& funcs) {
      std::string own_desc = "returns the description of a function, takes a string as an argument, for empty string returns the list of functions";
      auto descs = std::make_shared<std::unordered_map<std::string, std::string>>();
      (*descs)["help"] = own_desc;
      for (int i = 0; i < funcs.size(); i++) {
        (*descs)[funcs[i].name] = funcs[i].description;
      }
      return FunctionEntry(
        "help",
        [descs](std::shared_ptr<Thunk> s) -> GenericDataType {
            auto s_val = s -> force();
            if (s_val.myType != GenericDataType::STRING) {
              throw std::runtime_error("help takes a string as argument");
            }
            if (s_val.str_val == "") {
              std::string res = "";
              for (auto it = (*descs).begin(); it != (*descs).end(); ++it) {
                res += "\n" + it -> first;
              }
              return GenericDataType(res);
            }
            std::string res = (*descs)[s_val.str_val];
            return GenericDataType(res);
          },
        own_desc
        );
    }

};

class Parser {
  public:
    using FuncType = std::function<GenericDataType(std::shared_ptr<Thunk>)>;
    using FuncTypeAppl = std::pair<std::shared_ptr<Thunk>, std::shared_ptr<Thunk>>;

    class ParseReturnStruct {
      public:
        std::shared_ptr<GenericDataType> dataType = std::make_shared<GenericDataType>();
        std::vector<std::string>::iterator tokensStart;
        // std::vector<std::string>::iterator tokensEnd;

        ParseReturnStruct(
          GenericDataType& dT,
          std::vector<std::string>::iterator tS //,
          //std::vector<std::string>::iterator tE
          ) {
            dataType = std::make_shared<GenericDataType>(dT);
            tokensStart = tS;
            // tokensEnd = tE;
          }
        };

    // SYNTAX:
    // parentheses: (), [], {}
    // expr: funcdef, funcappl, variable, const
    // funcdef: lambda paramname . expr (abstraction)
    // funcappl: func param (application)
    // const: literal or const func
    // variable
    // M N P -> (M N) P
    // lambda x . M N -> lambda x . (M N)
    // lambda x y z . N -> lambda x . lambda y . lambda z . N
    // literals:
    // NIL: nil
    // BOOL: false, true
    // DOUBLE: same as normal, but accepts , instead of .
    // HASHTABLE: table - constructs an empty table

    static GenericDataType GenericDataTypeFromString(std::string s, std::vector<FunctionEntry> funcs) {
      std::vector<std::string> tokens = splitString(s, " \t\n", "()[]{}."); // split over whitespace symbols, also isolate the single char tokens

      // DEBUG OUTPUT OF TOKENS
      // std::cout << s << std::endl;
      /*
      for (int i = 0; i < tokens.size(); i++) {
        std::cout << i << " " << tokens[i] << "; ";
      }

      std::cout << std::endl;
      */
      ParseReturnStruct tmp = Parse(tokens.begin(), tokens.end(), {}, funcs);
      if (tmp.tokensStart != tokens.end()) {
        std::string unparsed = "";
        for (std::vector<std::string>::iterator it = tmp.tokensStart; it != tokens.end(); it++) {
          unparsed += *it + "; ";
        }
        throw std::runtime_error("some string remains unparsed: " + unparsed);
      }
      return *tmp.dataType;
    }

    static std::shared_ptr<GenericDataType> FunctionSequenceHelper(std::vector<std::shared_ptr<GenericDataType>> funVals) {
      if (funVals.size() == 0) {
        throw std::runtime_error("a sequence of 0 functions encountered");
      }
      // so, (((M N) P) Q) R
      if (funVals.size() == 1) {
        return funVals[0];
      }
      auto last = funVals[funVals.size() - 1];
      funVals.pop_back();
      auto first = FunctionSequenceHelper(funVals);
      if (first -> myType == GenericDataType::FUNCTION || first -> myType == GenericDataType::FUNCTIONAPPL || first -> myType == GenericDataType::VARIABLE) {
        GenericDataType tmp = GenericDataType();
        tmp.myType = GenericDataType::FUNCTIONAPPL;
        FuncTypeAppl res = FuncTypeAppl(std::make_shared<Thunk>(*first), std::make_shared<Thunk>(*last));
        tmp.funcAppl_val = std::make_shared<FuncTypeAppl>(res); // apply the function to the argument
        // throw std::runtime_error("function call not implemented yet");
        return std::make_shared<GenericDataType>(tmp);
      }
      throw std::runtime_error("tried to call a non-function type");
      // return first -> func_val(last); // simply call the function, what could possibly go wrong with unparsed varnames
    }

    static FuncType LambdaAbstractionHelper(std::string varname, std::shared_ptr<Thunk> expr) {
      // bind the variable
      // std::cout << "doing the bind of " + varname + " to " + expr.displayStr() << std::endl;
      if (expr->expr.myType == GenericDataType::HASHMAP ||
          expr->expr.myType == GenericDataType::NUMBER ||
          expr->expr.myType == GenericDataType::STRING ||
          expr->expr.myType == GenericDataType::BOOL ||
          expr->expr.myType == GenericDataType::NIL) {
        // we just return the body constant
        return [expr](std::shared_ptr<Thunk> x) -> GenericDataType { return expr -> expr; };
      }
      if (expr->expr.myType == GenericDataType::VARIABLE) {
        if (varname == expr->expr.var_val) {
          // we must return lambda x.x
          return [](std::shared_ptr<Thunk> x) -> GenericDataType { return x -> expr; };
        }
        // we just return the body constant
        return [expr](std::shared_ptr<Thunk> x) -> GenericDataType { return expr -> expr; };
      }
      if (expr->expr.myType == GenericDataType::FUNCTIONAPPL) {
        // we have lambda x.(f g)
        // we compute F = lambda x.f
        // we compute G = lambda x.G
        std::string var = varname;
        std::shared_ptr<Thunk> f_val = (expr->expr.funcAppl_val)->first;
        std::shared_ptr<Thunk> s_val = (expr->expr.funcAppl_val)->second;
        auto first = std::make_shared<FuncType>(LambdaAbstractionHelper(var, f_val));
        auto second = std::make_shared<FuncType>(LambdaAbstractionHelper(var, s_val));
        // lambda x.(f g) = lambda x.(F(x), G(x))
        return [first, second](std::shared_ptr<Thunk> x) -> GenericDataType {
          // std::cout << "called function application subsitution of varname, subsituting for: " + x.displayStr() << std::endl;
          // std::cout << "res of 1st is: " + (*first)(x).displayStr() << std::endl;
          // std::cout << "res of 2nd is: " + (*second)(x).displayStr() << std::endl;
          GenericDataType tmp = GenericDataType();
          tmp.myType = GenericDataType::FUNCTIONAPPL;
          FuncTypeAppl res = std::pair<std::shared_ptr<Thunk> , std::shared_ptr<Thunk> >(
            std::make_shared<Thunk>((*first)(x)),
            std::make_shared<Thunk>((*second)(x)));
          // std::cout << "success" << std::endl;
          tmp.funcAppl_val = std::make_shared<FuncTypeAppl>(res); // apply the function to the argument
          return tmp;
        };
      }

      if (expr->expr.myType == GenericDataType::FUNCTION) {
        // we have lambda x.f
        // f = lambda y.g
        // we return lambda x.(lambda y.LambdaAbstractionHelper(x, g(y)))
        // so
        // we return lambda x.(lambda y.LambdaAbstractionHelper(x, f(y)))
        std::shared_ptr<FuncType> f = expr->expr.func_val;
        std::string var = varname;
        std::string f_com = expr->expr.func_comment;
        return [var, f, f_com](std::shared_ptr<Thunk> x) -> GenericDataType {
          auto x_val = x;
          auto fun = [var, f, x_val](std::shared_ptr<Thunk> y) -> GenericDataType {
            GenericDataType tmp = (*f)(y);
            return LambdaAbstractionHelper(var, std::make_shared<Thunk>(tmp))(x_val);
          };
          GenericDataType tmp = GenericDataType();
          tmp.myType = GenericDataType::FUNCTION;
          //tmp.func_comment = "[" + var + "=" + x -> expr.displayStr() + "]" + f_com;
          size_t threshold = 50;
          std::string res_com = f_com.substr(0, std::min(f_com.size(), (size_t)threshold));
          if (f_com.size() > threshold) {
            res_com += "...";
          }
          tmp.func_comment = "[" + var + "]" + res_com;
          tmp.func_val = std::make_shared<FuncType>(fun);
          return tmp;
        };
      }
      throw std::runtime_error("unknown type in lambda abstraction");

    }

    static ParseReturnStruct LambdaSequenceHelper(std::vector<std::string> varnames, ParseReturnStruct rest) {
      if (varnames.size() == 0) {
        throw std::runtime_error("a sequence of 0 arguments in a lambda function declaration is invalid");
      }
      if (varnames.size() == 1) {
        GenericDataType tmp = GenericDataType();
        tmp.myType = GenericDataType::FUNCTION;
        tmp.func_val = std::make_shared<FuncType>(LambdaAbstractionHelper(varnames[0], std::make_shared<Thunk>(*rest.dataType)));
        tmp.func_comment = "lambda " + varnames[0] + "." + (*rest.dataType).displayStr();
        rest.dataType = std::make_shared<GenericDataType>(tmp);
        return rest;
      }
      std::string last = varnames[varnames.size() - 1];
      varnames.pop_back();
      return LambdaSequenceHelper(varnames, LambdaSequenceHelper({ last }, rest)); // do the last one, then pass back into self
    }


    static ParseReturnStruct Parse(std::vector<std::string>::iterator tokensStart, std::vector<std::string>::iterator tokensEnd, std::vector<std::string> validVarNames, std::vector<FunctionEntry>& functions) {
      if (tokensStart == tokensEnd) {
        throw std::runtime_error("tried to parse empty expression");
        GenericDataType tmp = GenericDataType();
        return ParseReturnStruct(tmp, tokensStart); // return nothing
      }

      // first check for parentheses

      // if we have closure parentheses, we must exit early
      std::string openingParentheses = "([{";
      std::string closureParentheses = ")]}";
      if (closureParentheses.find(*tokensStart) != std::string::npos) {
        GenericDataType tmp = GenericDataType();
        return ParseReturnStruct(tmp, tokensStart + 1); // return all remaining tokens
      }

      // check lambda
      if (*tokensStart == "lambda") {
        // lambda varname1 varname2 varname3 ... . exprUsingNewVarNames
        // ->
        // lambda varname1 . lambda varname2 . lambda varname3 ... . exprUsingVarNames
        // allow user to make nonsensical varnames for convenience, so just parse until dot
        std::vector<std::string> varnames;
        std::vector<std::string>::iterator tokenCurVarNames = tokensStart + 1;
        while (tokenCurVarNames != tokensEnd) {
          if (*tokenCurVarNames == ".") {
              break;
          }
          else {
            varnames.push_back(*tokenCurVarNames);
            validVarNames.push_back(*tokenCurVarNames);
            tokenCurVarNames++;
          }
        }
        if (tokenCurVarNames == tokensEnd) {
          throw std::runtime_error("lambda function arguments never ended with a '.'");
        }
        ParseReturnStruct tmp = Parse(tokenCurVarNames + 1, tokensEnd, validVarNames, functions);
        return LambdaSequenceHelper(varnames, tmp);
      }

      // check application
      // M -> M
      // M N -> M N
      // M N P -> (M N) P
      // so, try parse M N P Q R ... until lambda or ([{}])

      std::vector<std::shared_ptr<GenericDataType>> funVals;
      std::vector<std::string>::iterator tokenCur = tokensStart;
      while (tokenCur != tokensEnd) {
        int openType = 0;
        if (*tokenCur == "lambda") {
          // continue until lambda expr end
          ParseReturnStruct tmp = Parse(tokenCur, tokensEnd, validVarNames, functions);
          funVals.push_back(tmp.dataType);
          tokenCur = tmp.tokensStart;
        }
        else if ((openType = openingParentheses.find(*tokenCur)) != std::string::npos) {
          // continue until parenthesis end
          // if we have opening parentheses, we must evaluate the expression inside and check the end parentheses type
          ParseReturnStruct tmp = Parse(tokenCur + 1, tokensEnd, validVarNames, functions);
          if (closureParentheses.find(*tmp.tokensStart) != openType) {
            throw std::runtime_error("closing parentheses '" + *tmp.tokensStart + "' do not match opening parentheses '" + openingParentheses[openType] + "'");
          }
          ParseReturnStruct tmp2 = ParseReturnStruct(*tmp.dataType, tmp.tokensStart + 1); // return all remaining tokens, plus what we parsed inside the parentheses
          funVals.push_back(tmp2.dataType);
          tokenCur = tmp2.tokensStart;
        }
        else if (closureParentheses.find(*tokenCur) != std::string::npos) { break; }
        else {
          // add literal or variable or const or...
          // first check if we are a valid variable
          bool validVarName = false;
          for (int i = 0; i < validVarNames.size(); i++) {
            // order of check doesn't matter, so we can forward iterate
            if (validVarNames[i] == *tokenCur) {
              validVarName = true;
              break;
            }
          }
          if (validVarName) {
            // we are a variable name, simple as
            GenericDataType tmp = GenericDataType();
            tmp.myType = GenericDataType::VARIABLE;
            tmp.var_val = *tokenCur;
            funVals.push_back(std::make_shared<GenericDataType>(tmp));
          } else {
            // we are a literal or a builtin function, all builtin functions are to be added here

            int validFunName = -1;
            for (int i = 0; i < functions.size(); i++) {
              // order of check doesn't matter, so we can forward iterate
              if (functions[i].name == *tokenCur) {
                validFunName = i;
                break;
              }
            }
            if (validFunName != -1) {
                GenericDataType tmp = GenericDataType();
                tmp.myType = GenericDataType::FUNCTION;
                tmp.func_val = functions[validFunName].func_val;
                tmp.func_comment = functions[validFunName].name;
            funVals.push_back(std::make_shared<GenericDataType>(tmp));
            } else {
              // we are not a function, so try a literal
              if (*tokenCur == "nil") { // NIL
                funVals.push_back(std::make_shared<GenericDataType>());
              } else if (*tokenCur == "false") { // BOOL
                GenericDataType tmp = GenericDataType();
                tmp.myType = GenericDataType::BOOL;
                tmp.bool_val = false;
                funVals.push_back(std::make_shared<GenericDataType>(tmp));
              } else if (*tokenCur == "true") { // BOOL
                GenericDataType tmp = GenericDataType();
                tmp.myType = GenericDataType::BOOL;
                tmp.bool_val = true;
                funVals.push_back(std::make_shared<GenericDataType>(tmp));
              } else if ((*tokenCur)[0] == '\'') { // STRING
                if ((*tokenCur)[(*tokenCur).length() - 1] != '\'' || (*tokenCur).length() < 2) {
                  throw std::runtime_error("invalid string: " + *tokenCur);
                }
                GenericDataType tmp = GenericDataType((*tokenCur).substr(1, (*tokenCur).length() - 2));
                funVals.push_back(std::make_shared<GenericDataType>(tmp));
              } else if (*tokenCur == "table") { // HASHMAP
                GenericDataType tmp = GenericDataType();
                tmp.myType = GenericDataType::HASHMAP;
                tmp.hashmap_val = std::make_shared<std::unordered_map<GenericDataType, GenericDataType, GenericDataType::GenericDataTypeHash>>();
                funVals.push_back(std::make_shared<GenericDataType>(tmp));
              } else { // NUMBER
                std::string str = (*tokenCur);
                std::replace(str.begin(), str.end(), ',', '.');
                GenericDataType tmp = GenericDataType(std::stod(str));
                funVals.push_back(std::make_shared<GenericDataType>(tmp));
                // throw std::runtime_error("unexpected literal: " + *tokenCur);
              }
            }
          }
          tokenCur++;
        }
      }

      // we now have our function vals
      return ParseReturnStruct(*FunctionSequenceHelper(funVals), tokenCur);

      // throw std::runtime_error("unexpected token: " + *tokensStart);
      // return nothing
    }


};

class XMLLambdaParserInstance {
  public:
    using FuncTypeAppl = std::pair<std::shared_ptr<Thunk>, std::shared_ptr<Thunk>>;
    GenericDataType global;
    std::vector<FunctionEntry> funcs;
    bool debug_log_eval = false;

    XMLLambdaParserInstance() {
        // assign an empty hashmap to global
        global = GenericDataType();
        global.myType = GenericDataType::HASHMAP;
        global.hashmap_val =  std::make_shared<std::unordered_map<GenericDataType, GenericDataType, GenericDataType::GenericDataTypeHash>>();
        // set functions
        funcs = FunctionEntry::getDefaultFunctions();
    }

    // takes a few functions of global as input, and inserts the result into global at the result of name
    void AssignGlobal(std::string name, std::string expr) {
      GenericDataType parsedName = Parser::GenericDataTypeFromString(name, funcs);
      GenericDataType parsedExpr = Parser::GenericDataTypeFromString(expr, funcs);

      // run the functions with global as input
      GenericDataType nameRes = GenericDataType();
      nameRes.myType = GenericDataType::FUNCTIONAPPL;
      nameRes.funcAppl_val = std::make_shared<FuncTypeAppl>(std::make_shared<Thunk>(parsedName), std::make_shared<Thunk>(global));

      GenericDataType exprRes = GenericDataType();
      exprRes.myType = GenericDataType::FUNCTIONAPPL;
      exprRes.funcAppl_val = std::make_shared<FuncTypeAppl>(std::make_shared<Thunk>(parsedExpr), std::make_shared<Thunk>(global));

      // insert quick{
      (*global.hashmap_val)[nameRes.eval(debug_log_eval)] = exprRes.eval(debug_log_eval);
    }

    // takes a function of global and outputs the compiled variant, may be a function that can be used later
    GenericDataType Expression(std::string expr) {
      GenericDataType parsedExpr = Parser::GenericDataTypeFromString(expr, funcs);

      // run the function with global as input
      GenericDataType exprRes = GenericDataType();
      exprRes.myType = GenericDataType::FUNCTIONAPPL;
      exprRes.funcAppl_val = std::make_shared<FuncTypeAppl>(std::make_shared<Thunk>(parsedExpr), std::make_shared<Thunk>(global));

      return exprRes.eval(debug_log_eval);
    }
};