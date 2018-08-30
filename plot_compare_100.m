%plot-compare!

a = 1;

newParam = csvread('NewParam.csv');
fprintf('the new parameterisation is in green\n');
oldParam = csvread('OldParam.csv');
fprintf('the old parameterisation is in red\n');
tParam = csvread('tParam.csv');
fprintf('the t parameterisation is in blue\n');

NewScore = zeros(size(newParam(:,1),1),1);
OldScore = zeros(size(newParam(:,1),1),1);
tScore = zeros(size(newParam(:,1),1),1);



for i = 1:99
    
    New = ((newParam((i-1)*6+1,2:22)-newParam((i-1)*6+3,2:22)));
    Old= ((oldParam((i-1)*6+1,2:22)-oldParam((i-1)*6+3,2:22)));
    t = ((tParam((i-1)*6+1,2:22)-tParam((i-1)*6+3,2:22)));
    
    NewUpBig = New(New>1).^4;
    NewUp = (2*New(New<1&New>0)).^2;
    NewScoreDown = New(New<0).^2;
    NewScore(i) = sum(NewUpBig) + sum(NewUp) + sum(NewScoreDown);
    
    OldUpBig = Old(Old>1).^4;
    OldUp = (2*New(Old<1&Old>0)).^2;
    OldScoreDown = Old(Old<0).^2;
    OldScore(i) = sum(OldUpBig) + sum(OldUp) + sum(OldScoreDown);    

    tUpBig = t(t>1).^4;
    tUp = (2*t(t<1&t>0)).^2;
    tScoreDown = t(t<0).^2;
    tScore(i) = sum(tUpBig) + sum(tUp) + sum(tScoreDown);    
end

hold all
plot(1:99,NewScore(1:99),'r')
plot(1:99,OldScore(1:99),'g')
%plot(1:99,tScore(1:99),'b')
hold off

%%%%

while 1==0%a < 100
    if a == 38
        fprintf('graph 38 does not exist for some reason\n')
        fprintf('displaying graph 39...\n')
    end
    if a > 38
        a=a-1;
    end
    close all
    hold all
    new1=plot(newParam((a-1)*6+1,2:22),'--g');
    ylabel('Generation, kWH');
    title(['House number ',num2str(a)])
    new2=plot(newParam((a-1)*6+2,2:22),':g');
    new2.LineWidth=2;
    new3=plot(newParam((a-1)*6+3,2:22),'g');
    new4=plot(newParam((a-1)*6+4,2:22),'g');
    new5=plot(newParam((a-1)*6+5,2:22),'g');
    new6=plot(newParam((a-1)*6+6,2:22),'g');

    old1=plot(oldParam((a-1)*6+1,2:22),'--r');
    old2=plot(oldParam((a-1)*6+2,2:22),':r');
    old2.LineWidth=2;
    old3=plot(oldParam((a-1)*6+3,2:22),'r');
    old4=plot(oldParam((a-1)*6+4,2:22),'r');
    old5=plot(oldParam((a-1)*6+5,2:22),'r');
    old6=plot(oldParam((a-1)*6+6,2:22),'r');
    
    t1=plot(tParam((a-1)*6+1,2:22),'--b');
    t2=plot(tParam((a-1)*6+2,2:22),':b');
    t2.LineWidth=2;
    t3=plot(tParam((a-1)*6+3,2:22),'b');
    t4=plot(tParam((a-1)*6+4,2:22),'b');
    t5=plot(tParam((a-1)*6+5,2:22),'b');
    t6=plot(tParam((a-1)*6+6,2:22),'b');
    hold off
    a = input('what graph would you like to see next? ');
    
end