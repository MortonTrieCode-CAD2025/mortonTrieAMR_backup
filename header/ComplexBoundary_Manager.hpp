/**
 * @file SolidStream_Manager.h
 * @brief fluid stream aroung obstacle
 * @version 0.1
 * @date 2024-02-26
 * 
 */

#pragma once

#include "user.h"

class ComplexBoundary_Model
{
protected:
    User* user_;

public:
    ComplexBoundary_Model(User* user) : user_(user) {}

    virtual void treat() = 0;

    virtual ~ComplexBoundary_Model() = default;
};

class UniformBoundary : public ComplexBoundary_Model
{
public:
    UniformBoundary(User* user) : ComplexBoundary_Model(user) {}

    void treat() override;
};

class SinglePointInterpolation : public ComplexBoundary_Model
{
public:
    SinglePointInterpolation(User* user) : ComplexBoundary_Model(user) {}

    void treat() override;
};