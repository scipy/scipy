obj1 = simpleclass;
obj2 = simpleclass;
obj3 = anotherclass;
obj4 = anotherclass;
obj1.any_field_3 = obj1;
obj2.char_field_1 = 'foo';
obj2.any_field_3 = [obj3, obj4];
save -v6 testclass_v6 obj1 obj2 obj3;
save -v7 testclass_v7 obj1 obj2 obj3;
