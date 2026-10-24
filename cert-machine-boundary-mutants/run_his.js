// run cert-machine's own circles decider (instruments/easota/decide.js, the ledger's path) on each mutant
'use strict';
const fs = require('fs'), path = require('path');
const ROOT = process.argv[2], DIR = process.argv[3];
const L = require(path.join(ROOT, 'instruments/easota/lib.js'));
const D = require(path.join(ROOT, 'instruments/easota/decide.js'));
const out = {};
// usage: node run_his.js <cert-machine root> <this directory>
for (const f of fs.readdirSync(DIR).filter((n) => /^circles__.*\.json$/.test(n))) {
  const [row, id] = [f.split('__').slice(0, 2).join('__'), f.split('__')[2].replace('.json', '')];
  const C = JSON.parse(fs.readFileSync(path.join(DIR, f), 'utf8')).circles.map((c) => c.map(L.fromJsonNumber));
  const r = D.circles(C);
  out[row + '/' + id] = { word: r.witnessed ? 'WITNESSED' : (r.repair ? 'REPAIRED' : 'UNWITNESSED'), overlaps: r.overlaps, boxSlackSign: L.Q.sign(r.boxSlack), worst: [r.worstPair.i, r.worstPair.j] };
}
console.log(JSON.stringify(out, null, 1));
