// public/js/pubchem.js with a stand-in for fetch: the page never needs the network in tests.
import assert from "node:assert/strict";
import { test } from "node:test";

import { PubChemError, countRows, fetchSearch } from "../js/pubchem.js";

const answer = (status, body, statusText = "") => async (url) => ({
  ok: status >= 200 && status < 300,
  status,
  statusText,
  url,
  text: async () => body,
});

test("fetchSearch returns PubChem's text", async () => {
  const calls = [];
  const fetchFn = async (url) => {
    calls.push(url);
    return answer(200, '"CID","ConnectivitySMILES"\n')(url);
  };
  assert.equal(await fetchSearch("https://pubchem.example/x", fetchFn), '"CID","ConnectivitySMILES"\n');
  assert.deepEqual(calls, ["https://pubchem.example/x"]);
});

test("fetchSearch reports PubChem's status when it refuses", async () => {
  await assert.rejects(
    fetchSearch("u", answer(503, "busy", "Service Unavailable")),
    (error) => error instanceof PubChemError && /503/.test(error.message) && /Service Unavailable/.test(error.message),
  );
});

test("fetchSearch reports an unreachable PubChem", async () => {
  const offline = async () => {
    throw new TypeError("Failed to fetch");
  };
  await assert.rejects(
    fetchSearch("u", offline),
    (error) => error instanceof PubChemError && /could not reach PubChem/.test(error.message) && /Failed to fetch/.test(error.message),
  );
});

test("countRows counts data lines after the header", () => {
  assert.equal(countRows('"CID","ConnectivitySMILES"\n1,"a"\n2,"b"\n'), 2);
  assert.equal(countRows('"CID","ConnectivitySMILES"\r\n1,"a"\r\n'), 1);
  assert.equal(countRows('"CID","ConnectivitySMILES"\n'), 0);
  assert.equal(countRows(""), 0);
  assert.equal(countRows('"CID","ConnectivitySMILES"\n1,"a"\n\n'), 1);
});
