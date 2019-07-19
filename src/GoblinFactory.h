#ifndef GOBLIN_FACTORY_H
#define GOBLIN_FACTORY_H

#include <map>
#include <string>
#include "GoblinUtils.h"
#include "GoblinParamSet.h"

namespace Goblin {

    template<typename Product, typename P1 = NullType,
        typename P2 = NullType, typename P3 = NullType>
    class Creator {
    public:
        virtual ~Creator() {}
        virtual Product* create() const { return NULL; }
        virtual Product* create(P1 p1) const { return NULL; };
        virtual Product* create(P1 p1, P2 p2) const { return NULL; };
        virtual Product* create(P1 p1, P2 p2, P3 p3) const { return NULL; };
    };


    template<typename Product, typename P1 = NullType, 
        typename P2 = NullType, typename P3 = NullType>
    class Factory {
    public:
        ~Factory();
        bool registerCreator(const string& type, 
            Creator<Product, P1, P2, P3>* creator);
        Product* create(const string& type);
        Product* create(const string& type, P1 p1);
        Product* create(const string& type, P1 p1, P2 p2);
        Product* create(const string& type, P1 p1, P2 p2, P3 p3);
        bool setDefault(const string& type);

    private:
        typedef map<string, Creator<Product, P1, P2, P3>* > CreatorTable;
        CreatorTable mCreatorTable;
        string mDefaultType;
    };


    template<typename Product, typename P1, typename P2, typename P3> 
    Factory<Product, P1, P2, P3>::~Factory() {
        typename CreatorTable::iterator it = mCreatorTable.begin();
        for (; it != mCreatorTable.end(); ++it) {
            delete it->second;
        }
        mCreatorTable.clear();
    }

    template<typename Product, typename P1, typename P2, typename P3> 
    bool Factory<Product, P1, P2, P3>::registerCreator(const string& type, 
        Creator<Product, P1, P2, P3>* creator) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            cerr << "type " << type << " already registered" << endl;
            return false;
        } else {
            mCreatorTable[type] = creator;
            return true; 
        }
    }


    template<typename Product, typename P1, typename P2, typename P3> 
    Product* Factory<Product, P1, P2, P3>::create(const string& type) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            return mCreatorTable[type]->create(); 
        } else if (mCreatorTable.find(mDefaultType) != mCreatorTable.end()) {
            cerr << "unrecognized type: " << type << 
                ", fallback to " << mDefaultType <<endl;
            return mCreatorTable[mDefaultType]->create();
        } else {
            cerr << "unrecognized type: " << type << endl;
            return NULL;
        }
    }

    template<typename Product, typename P1, typename P2, typename P3> 
    Product* Factory<Product, P1, P2, P3>::create(const string& type,
        P1 p1) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            return mCreatorTable[type]->create(p1); 
        } else if (mCreatorTable.find(mDefaultType) != mCreatorTable.end()) {
            cerr << "unrecognized type: " << type << 
                ", fallback to " << mDefaultType <<endl;
            return mCreatorTable[mDefaultType]->create(p1);
        } else {
            cerr << "unrecognized type: " << type << endl;
            return NULL;
        }
    }

    template<typename Product, typename P1, typename P2, typename P3> 
    Product* Factory<Product, P1, P2, P3>::create(const string& type,
        P1 p1, P2 p2) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            return mCreatorTable[type]->create(p1, p2); 
        } else if (mCreatorTable.find(mDefaultType) != mCreatorTable.end()) {
            cerr << "unrecognized type: " << type << 
                ", fallback to " << mDefaultType <<endl;
            return mCreatorTable[mDefaultType]->create(p1, p2);
        } else {
            cerr << "unrecognized type: " << type << endl;
            return NULL;
        }
    }

    template<typename Product, typename P1, typename P2, typename P3> 
    Product* Factory<Product, P1, P2, P3>::create(const string& type,
        P1 p1, P2 p2, P3 p3) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            return mCreatorTable[type]->create(p1, p2, p3); 
        } else if (mCreatorTable.find(mDefaultType) != mCreatorTable.end()) {
            cerr << "unrecognized type: " << type << 
                ", fallback to " << mDefaultType <<endl;
            return mCreatorTable[mDefaultType]->create(p1, p2, p3);
        } else {
            cerr << "unrecognized type: " << type << endl;
            return NULL;
        }
    }

    template<typename Product, typename P1, typename P2, typename P3> 
    bool Factory<Product, P1, P2, P3>::setDefault(const string & type) {
        if (mCreatorTable.find(type) != mCreatorTable.end()) {
            mDefaultType = type;
            return true;
        } else {
            cerr << "fail to set default creator to " << type <<
                ": creator not registered " << endl;
            return false;
        }
    }
}

#endif //GOBLIN_FACTORY_H
