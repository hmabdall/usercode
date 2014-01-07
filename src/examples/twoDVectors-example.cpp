#include <iostream>
#include <vector>
using namespace std;

int main()
{   
   int numberOfEntries = 3;
   
    vector<vector<double> > event;
    
    for (int i = 0; i < numberOfEntries - 1; i++){
    
        event.push_back(vector<double>());
        
        
    }
    
       
           cout << "Hello World bitchez" << endl; 

    for (int jentry = 0; jentry < numberOfEntries - 1; jentry++)
    
    {
        event[jentry].push_back(28.5);
        event[jentry].push_back(78.5);
    }    
    

    for (int jentry = 0; jentry < numberOfEntries - 1; jentry++){
        cout << "event no:" << jentry << " " << event[jentry][0] << endl;
        cout << "event no:" << jentry << " "<< event[jentry][1] << endl;
    }
    
        
    
    
    
   
   return 0;
}


