def c
pipeline {
  agent any

  stages {
    stage('Build') {
      steps {
        script {
          c = docker.build("phewas-development/betamatch:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
          docker.withRegistry('http://gcr.io/phewas-development', 'gcr:phewas-development') {
            c.push("test-${env.BUILD_NUMBER}")
          }
    }
      }
    }
    stage('Tests') {
      /*set up tests*/
      steps{
        script{
          c.inside("-u root"){sh """python3 -m pip install pytest
                                    cd /usr/local/betamatch
                                    python3 -m pytest
                                    wget https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar
                                    java -jar womtool-48.jar wdl/betamatch_github.wdl -i wdl/betamatch_github.wdl"""}
        }
        /*sh 'python --version'
        sh 'python3 -m pip install pylint pytest safety pyflakes mypy prospector bandit'
        sh 'curl -O https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar'
        sh 'chmod +x womtool-48.jar'
        sh './womtool-48.jar betamatch/wdl/betamatch_github.wdl -i betamatch/wdl/betamatch_github.wdl'
        sh 'cd betamatch'
        sh 'python3 -m pytest'*/
      }
    }
    stage('Metrics') {
      steps {
        withSonarQubeEnv('sonar') {
          sh "${tool("sonar")}/bin/sonar-scanner \
          -Dsonar.projectKey=${JOB_NAME} \
          -Dsonar.sources=. \
          -Dsonar.css.node=. \
          -Dsonar.host.url=${DEFAULT_SONAR_URL} \
          -Dsonar.login=${SONAR_LOGIN}"
        }
      }
    }
  }
}