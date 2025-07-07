// /**
//  * @file ZouHe_Pressure.h
//  * @brief Non-equilibrum Bounce-Back (ZouHe) Pressure BC = Neumann BC (suitable for Low Re Number)
//  * @date 2023-07-18
//  */

// #pragma once
// #include "user.h"
// #include "D3Q19_BC_Manager.h"

// class nonEquilibrumBounceBack_Pressure_BC : public D3Q19_BC_Strategy
// {
// public:
//     virtual void applyBCStrategy() = 0;
//     virtual void initialBCStrategy() = 0;
// };

// class nonEquilibrumBounceBack_Pressure_West : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_West() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };

// class nonEquilibrumBounceBack_Pressure_East : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_East() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };

// class nonEquilibrumBounceBack_Pressure_North : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_North() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };

// class nonEquilibrumBounceBack_Pressure_South : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_South() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };

// class nonEquilibrumBounceBack_Pressure_Top : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_Top() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };

// class nonEquilibrumBounceBack_Pressure_Bot : public nonEquilibrumBounceBack_Pressure_BC
// {
// public:
//     nonEquilibrumBounceBack_Pressure_Bot() = default;
//     void applyBCStrategy() override;
//     void initialBCStrategy() override;

// private:
//     /**
//      * @todo double pressure/density
//      */
// };