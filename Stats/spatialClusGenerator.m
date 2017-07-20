function [clusMass, clusterIdx] = spatialClusGenerator(data,perm,alpha,elecDists_idx)


clusMass = 0;
for i = 1:size(data,1) % iterate through first dimension should be channels
    for j = 1:size(data,2) % iterate through second dimension should be time
        if size(data,2)>1 % there is a time component
            foo = squeeze(perm(i,j,:));
            p_val(i,j) = sum(foo>data(i,j))/length(foo);
        else % there is no time component
            foo = perm(i,:);
            p_val(i,j) = sum(foo>data(i))/length(foo);
        end
    end
end

onoff = p_val2<=alpha;

onoff = permute(onoff,[1,3,2]); % format: onoff(Chan,Freq,Time)
[clusPos,num] = findcluster(onoff,elecDists_idx);
clusPos = squeeze(clusPos);

for j = 1:num
    idx = find(clusPos==j);
    clusMass(j) =  sum(data(idx));
%     clusMass(j) =  length(idx);
    [I,J] = ind2sub(size(clusPos),idx);
    clusterIdx{j} = [I,J];
end