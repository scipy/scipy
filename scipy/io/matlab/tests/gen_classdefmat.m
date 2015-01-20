obj1 = simpleclass;
obj2 = simpleclass;
obj3 = anotherclass;
obj2.char_field_1 = 'foo';
obj2.any_field_3 = obj3;
save -v6 testclass obj1 obj2 obj3;
