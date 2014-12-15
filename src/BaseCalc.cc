#include "LJMet/Com/interface/BaseCalc.h"
#include "LJMet/Com/interface/LjmetEventContent.h"

BaseCalc::BaseCalc():
mName(""),
mLegend("")
{
}

void BaseCalc::SetValue(std::string name, bool value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, int value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, double value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, std::vector<bool> value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, std::vector<int> value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::SetValue(std::string name, std::vector<double> value)
{
    std::string _name = name + "_" + mName;
    mpEc->SetValue(_name, value);
}

void BaseCalc::init()
{
    mLegend = "[" + mName + "]: ";
    std::cout << mLegend << "registering " << mName << std::endl;
}
