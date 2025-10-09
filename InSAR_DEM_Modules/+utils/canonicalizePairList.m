function canon = canonicalizePairList(pairList)
%CANONICALIZEPAIRLIST Ensures (ii1,jj1)-(ii2,jj2) are always ordered consistently
canon = pairList;
for i = 1:size(pairList,1)
    ii1 = pairList(i,1); jj1 = pairList(i,2);
    ii2 = pairList(i,3); jj2 = pairList(i,4);
    if ii1 > ii2 || (ii1 == ii2 && jj1 > jj2)
        canon(i,:) = [ii2 jj2 ii1 jj1];
    end
end
canon = unique(canon, 'rows');
end