function User_Active_Set = User_Active_Set_Generate(Q, TA, CA, TS)
%-----------------------------�����û���Ծλ��--------------------------  
%�����������ÿ��ʱ϶��Ծ�û�������TA (Total acttive)��
% ÿ��ʱ϶��ͬ�Ļ�Ծ�û�����CA (Common active)��TA-CA����Ծ�û��ǲ�ͬ�ġ�

DA = TA-CA; %��̬�仯�û�����(Dynamic active)
Set_Random=zeros(DA,TS); %���ö�̬�û�����
User_Active_Set=zeros(TA,TS);%���û�Ծ�û�����

User_Index = randperm(Q); 
Set_Common=User_Index(1:CA);%������Ծ�û�

for i_TS = 1 : TS
      random_index=randperm(Q-CA);  %ͬһʱ϶�û�����һ�����
   for i_DA = 1 : DA
      Set_Random(i_DA,i_TS)=User_Index(random_index(i_DA)+CA); %������̬��Ծ�û�����
   end
   User_Active_Set(:,i_TS)=[Set_Common';Set_Random(:,i_TS)]; %������Ծ�û����ϣ���������̬��  
end

end
