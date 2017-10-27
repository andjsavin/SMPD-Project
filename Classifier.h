#include "database.h"

class Classifier
{
	
	public:
		Classifier( const database &trainSet, const database &testSet ) : trainingSet(trainSet), testSet(testSet){}
		
		virtual void train() = 0;
		
		virtual float test() = 0; 
		
	private:
	
		database trainingSet;
		database testSet;
		
	
}