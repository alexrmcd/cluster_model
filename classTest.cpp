/////// 		classTest.cpp		4-21-16 ////////

#include <iostream>

class Hotel{
public:
	int rooms;
	int employees;
	int guests;
	double pproom;
	double wage ; //per day 
	double income ; //in millions


	Hotel() : rooms(120) , employees(100), guests(10) , pproom(150) , wage (100), income (0) {
		std::cout << "Constructing Hotel ...." << std::endl;
	}


	void newguests(int newguests){


		if (guests + newguests <= rooms)
			guests  = guests + newguests;
		
		else if (guests + newguests >= rooms)
			std::cout << "Fully booked!!" << std::endl;
					
	}


	int getguests(){
		return guests;
	}







};


main(){
	Hotel california;

	california.newguests(30);




	std::cout << california.getguests() << std::endl;

}