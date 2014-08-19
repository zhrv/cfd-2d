#pragma once
class LimiterDG
{
public:
	LimiterDG();
	virtual ~LimiterDG();
	virtual void run() = 0;
};

