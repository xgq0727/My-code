function User_Active_Set = User_Active_Set_Generate(Q, TA, CA, TS)
%-----------------------------设置用户活跃位置--------------------------  
%仿真中设的是每个时隙活跃用户数都是TA (Total acttive)，
% 每个时隙共同的活跃用户数是CA (Common active)，TA-CA个活跃用户是不同的。

DA = TA-CA; %动态变化用户数量(Dynamic active)
Set_Random=zeros(DA,TS); %设置动态用户集合
User_Active_Set=zeros(TA,TS);%设置活跃用户集合

User_Index = randperm(Q); 
Set_Common=User_Index(1:CA);%公共活跃用户

for i_TS = 1 : TS
      random_index=randperm(Q-CA);  %同一时隙用户集合一起产生
   for i_DA = 1 : DA
      Set_Random(i_DA,i_TS)=User_Index(random_index(i_DA)+CA); %产生动态活跃用户集合
   end
   User_Active_Set(:,i_TS)=[Set_Common';Set_Random(:,i_TS)]; %产生活跃用户集合（公共＋动态）  
end

end
