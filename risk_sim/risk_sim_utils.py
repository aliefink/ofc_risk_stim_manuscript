import numpy as np
import random
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.optimize import minimize # minimize function is used for parameter recovery 
import seaborn as sns 
from scipy.stats import pearsonr
import statsmodels.api as sm


def simulate_gamble_risk(risk_aversion,task_info,rep):
    #inputs: 
    #risk - risk aversion 
    #task_info - example task structure from s02 data
    #rep - number of times to run simulation
    #trials - number of trials for simulation (not included for gamble always 200)

    rep_list = []
    tr = []
    choice_prob = []
    choice_pred = []
    util_g = []
    util_s = []
    p_g = []
    p_s = []
    safe = []
    high = []
    low = []
    win_p = []
    payoffs = []
    outcome = []
    total_profit = []
    

    for rep in range(rep):
        profit = 0
        trials = task_info['round']-1
        #loop through trials
        for trial in trials:

            safe_bet = 10 #safebet always $10 in gambling task
            low_bet = 0 #lowbet (gamble loss) always $0 in gambling task
            high_bet = task_info['pot.payoff'].iloc[trial] #potential payoff for every trial is the highbet offer equivalent
            win_prob = task_info['win.prob'].iloc[trial] #win probability is different for every trial (=(1-(num_shown/100))


            safe.append(safe_bet)
            high.append(high_bet)
            low.append(low_bet)
            win_p.append(win_prob)

            # transform to high bet value to utility (gamble)
            weighted_high_bet = win_prob * ((high_bet)**risk_aversion)
            util_gamble = weighted_high_bet + low_bet #unnecessary step because low is always 0 but included for clarity
            # util_safe = safe_bet #this form assumes safe is not weighted by risk aversion because always guaranteed
            util_safe = (safe_bet)**risk_aversion #this is an option but probably not what ignacio wants

            #choice probability 
            inverse_temp = 5 #temp param is fixed to 1 (neither exploratory or exploitative) - not important for this simulation

            # convert EV to choice probabilities via softmax
            p_gamble = np.exp(inverse_temp*util_gamble) / ( np.exp(inverse_temp*util_gamble) + np.exp(inverse_temp*util_safe) )
            p_safe = np.exp(inverse_temp*util_safe) / ( np.exp(inverse_temp*util_gamble) + np.exp(inverse_temp*util_safe) )

            if np.isnan(p_gamble): #if utility is too large, cannot run prob calculation bc denom too large
                p_gamble = 0.99
                p_safe = 0.01 #avoid absolute decisions
            if np.isnan(p_safe):
                p_safe = 0.99
                p_gamble = 0.01 #avoid absolute decisions 
            

            util_g.append(util_gamble)
            util_s.append(util_safe)
            p_g.append(p_gamble)
            p_s.append(p_safe)

            #choice = random.choices(['gamble','safe'],weights=[p_gamble,p_safe])[0]
            if p_gamble > p_safe:
                choice = 'gamble'
            else:
                choice = 'safe'
            choice_pred.append(choice)

            if choice == 'gamble':
                choice_prob.append(p_gamble)
                gamble_options = [high_bet, low_bet]
                payoff = np.random.choice(gamble_options,1,p=[win_prob,1-win_prob])[0]
                payoffs.append(payoff)
                if payoff == low_bet:
                    outcome.append('gamble_loss')
                else:
                     outcome.append('gamble_win')

            else:
                choice_prob.append(p_safe)
                payoff = 10
                payoffs.append(payoff)
                outcome.append('safe')

            profit = profit + payoff 
            total_profit.append(profit)
            tr.append(trial)
            rep_list.append(rep) 


    data = {'rep':rep_list,'tr':tr,'ChoicePred':choice_pred,'ChoiceProb':choice_prob,'util_gamble':util_g,'util_safe':util_s,'p_gamble':p_g,'p_safe':p_s,
                       'SafeBet':safe,'HighBet':high,'LowBet':low,'win_prob':win_p,'payoff':payoffs,'outcome':outcome,'total_profit':total_profit}
    DF = pd.DataFrame(data)
    
    return DF






def simulation_norm_gamble_choices(df): #to-do input column names to make this robust to standard + util 

    dict = {}
    offer_norm = df['util_gamble']/df['util_safe']
    quant = np.quantile(offer_norm,q=(0,0.2,0.4,0.6,0.8,1),axis=0)
    x_axis = [np.mean(quant[i:i+2],dtype=np.float64) for i in range(5)]
    dec = np.array(df['ChoicePred'].replace(['gamble','safe'],[1,0]))
    dec[(dec != 1) & (dec !=2)] = 0
    gamble_zip = list(zip(offer_norm,dec))
    dict['norm_evs'] = np.array(offer_norm)
    dict['choices'] = np.array(dec)
    dict['x_axis'] = x_axis
    norm_range = []
    choice_props = []
    for r in range(len(quant)-1):
        ev_range = np.array([quant[r],quant[r+1]])
        gamble_count = [z[1] for z in gamble_zip if z[0] >= ev_range[0] and z[0] <= ev_range[1]]
        ev_num = sum(gamble_count)
        ev_prop = ev_num/len(gamble_count)
        norm_range.append(ev_range)
        choice_props.append(ev_prop)
    dict['norm_range'] = np.array(norm_range)
    dict['choice_props'] = np.array(choice_props)
    
    
    return dict