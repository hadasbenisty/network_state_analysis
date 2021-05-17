function t_win = get_sliding_window(winhopS, winsizeS, L, fs)


winst = 1:winhopS*fs:L;
winen = winst+winsizeS*fs;
winst=winst(winen<L);
winen=winen(winen<L);
t_win = (winst+winen)/2;

