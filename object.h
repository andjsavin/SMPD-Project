#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>

class Object
{
private:
	int classID;
    std::string className;
    std::vector<double> features;


public:

    Object(const std::string &className, const std::vector<double> &features) :classID(-1), className(className), features(features)
    {
    }

    std::string getClassName() const;
    size_t getFeaturesNumber() const;
    const std::vector<double> &getFeatures() const;
};



#endif // OBJECT_H
