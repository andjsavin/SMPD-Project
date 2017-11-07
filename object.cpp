#include "object.h"

    std::string Object::getClassName() const
    {
        return className;
    }

    size_t Object::getFeaturesNumber() const
    {
        return features.size();
    }

    const std::vector<double> &Object::getFeatures() const
    {
        return features;
    }



