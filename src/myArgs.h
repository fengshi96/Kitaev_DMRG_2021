//
// Created by shifeng on 7/30/21.
//

#ifndef ITENSOR_2D_MYARGS_H
#define ITENSOR_2D_MYARGS_H

#include "itensor/util/args.h"
#include "itensor/types.h"

class myArgs : public itensor::Args {
public:
    template<typename T>
    void add(Name const& name, std::vector<T>) {

    }

    template<typename T>
    std::vector<T> getVector(Name const& name) const {
        //if(defined(name)) return get(name).vectorVal();
        return {0};
    }

    template<typename T>
    class myVal
            {
            public:
                enum Type { Vector };
            private:
                Name name_;
                Type type_;
                std::string sval_;
                std::vector<T> vec_;
                itensor::Real rval_;
            public:


                myVal();
                myVal(const char* name);
                myVal(Name const& name);
                myVal(Name const& name, std::vector<T> const& vec_) {}

                std::vector<T> vectorVal() const {return vec_; }




            };
};



#endif //ITENSOR_2D_MYARGS_H
