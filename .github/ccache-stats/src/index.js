import * as core from "@actions/core";
import * as github from "@actions/github";

const { exec } = require('child_process');


try {
  // `who-to-greet` input defined in action metadata file
  const nameToGreet = core.getInput("who-to-greet");
  core.info(`Hello ${nameToGreet}!`);

  exec('ccache -s', (error, stdout, stderr) => {
    if (error) {
      core.setOutput("stats", `exec error: ${error.message}`);
      return;
    }
    if (stderr) {
      core.setOutput("stats", `stderr: ${stderr}`);
      return;
    }
    code.setOutput("stats", `stdout: ${stdout}`);
  });

  // Get the JSON webhook payload for the event that triggered the workflow
  //const payload = JSON.stringify(github.context.payload, undefined, 2);
  //core.info(`The event payload: ${payload}`);

} catch (error) {
  core.setFailed(error.message);
}
