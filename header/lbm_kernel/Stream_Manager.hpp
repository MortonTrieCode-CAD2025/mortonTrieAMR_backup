/**
 * @file Stream_Model.h
 * @brief Define collision model for the LBM kernel in the factory module
 * @version 0.1
 * @date 2024-02-25
 * 
 */

#pragma once

#include "user.h"

class Stream_Model
{
protected:
    User* user_;

public:
    Stream_Model(User* user) : user_(user) {}

    virtual void stream(D_int level) = 0;

    virtual ~Stream_Model() = default;
};

class ABPatternPush : public Stream_Model
{
public:
    ABPatternPush(User* user) : Stream_Model(user) {}

    void stream(D_int level) override;
};

class ABPatternPull : public Stream_Model
{
public:
    ABPatternPull(User* user) : Stream_Model(user) {}

    void stream(D_int level) override;
};

// class ABCollisionStreamFusionPush : public Stream_Model
// {
// public:
//     ABCollisionStreamFusionPush(User* user) : Stream_Model(user) {}

//     void stream(D_int level) override;
// };

// class ABCollisionStreamFusionPull : public Stream_Model
// {
// public:
//     ABCollisionStreamFusionPull(User* user) : Stream_Model(user) {}

//     void stream(D_int level) override;
// };

// class AAPattern : public Stream_Model
// {
// public:
//     AAPattern(User* user) : Stream_Model(user) {}
    
//     void stream(D_int level) override;
// };