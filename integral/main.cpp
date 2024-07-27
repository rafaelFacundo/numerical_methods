/* 



static vector<Matrice> rangeKuttaThirdOrderMethod_alternative(Matrice So, double deltaT, Matrice (*function)(Matrice), int numberOfStates)
    {
        vector<Matrice> states = {So};
        Matrice Si_bar, Si_bar_half, Si = Matrice();
        Matrice F1, F2, F3 = Matrice();
        for (int i = 1; i <= numberOfStates; ++i)
        {
            Si_bar_half = states[i-1] + function(states[i-1]) * (deltaT/2);
            Si_bar = states[i-1] + (function(Si_bar_half)*2 - function(states[i-1])) * deltaT;
            Si = states[i-1] + ( function(states[i-1]) * (1.0/6.0) + function(Si_bar_half) * (4.0/6.0) + function(Si_bar) * (1.0/6.0) ) * deltaT;
            states.push_back(Si);
        }
        return states;
    } 
    
    
    
    
*/