#include<iostream>
using namespace std ;


void Label_Maker(string * string1, int string_length)
{

	int i = 0;
	
	for (i = string_length; i >= (string_length - 4); i--)
	{
		string1->erase(i) ;
	}

	
} 

