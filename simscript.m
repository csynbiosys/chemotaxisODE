Lt=[0:2000];
Lv=100+500.*(Lt>500)+500.*(Lt>1000)+100000.*(Lt>1500);

[t,y] = ode23s(@(t,y) chemosig(t,y,Lt,Lv),[0 2000],[2.085]); assemblea; plot(t,y,t,a,Lt,(Lv-100)./max(Lv));