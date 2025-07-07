/**
 * @file ourflow.h
 * @brief 
 * @date 2023-07-24
 */
#define BC_NO_GHOST
#ifndef BC_NO_GHOST
#pragma once
#include "user.h"
#include "D3Q19_BC_Manager.hpp"

class outflow_BC : public D3Q19_BC_Strategy
{
public:
	virtual void applyBCStrategy() = 0;
	virtual void initialBCStrategy() = 0;
};

class outflow_West : public outflow_BC
{
public:
	outflow_West() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};

class outflow_East : public outflow_BC
{
public:
	outflow_East() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};

class outflow_South : public outflow_BC
{
public:
	outflow_South() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};

class outflow_North : public outflow_BC
{
public:
	outflow_North() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};

class outflow_Bot : public outflow_BC
{
public:
	outflow_Bot() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};

class outflow_Top : public outflow_BC
{
public:
	outflow_Top() = default;
	void applyBCStrategy() override;
	void initialBCStrategy() override;
};
#endif