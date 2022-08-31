
const fs = require('fs');

let txt = fs.readFileSync('sd22/data.json');
let s = JSON.parse(txt);
console.log(s);
