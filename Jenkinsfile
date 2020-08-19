//global define of container variable so that multiple stages can use the build container
def c
//pipeline start
pipeline {
  //agent, i.e. which machine or container can execute this pipeline.
  agent any
  //Block for different stages of the pipeline
  stages {
    //build stage. Build the container and send it to GCR
    stage('Build') {
      steps {
        //stage consists of a script...
        script {
          //building the container. Note that the dockerfile and context are specified in second argument
          //this is necessary so that dockerfile does not need to be in repo root
          c = docker.build("phewas-development/betamatch:test-" + "$BUILD_NUMBER", "-f docker/Dockerfile ./")
          //push docker to registry, only phewas-development is configured right now 
          docker.withRegistry('http://gcr.io/phewas-development', 'gcr:phewas-development') {
            c.push("test-${env.BUILD_NUMBER}")
          }
    }
      }
    }
    //testing stage, insert your tests here
    stage('Tests') {
      steps{
        script{
          //just run the tests inside the container built in last step.
          //This is currently the easiest solution to set up the testing environment that I know of
          c.inside("-u root"){sh """python3 -m pip install pytest
                                    cd /usr/local/betamatch
                                    python3 -m pytest
                                    wget https://github.com/broadinstitute/cromwell/releases/download/48/womtool-48.jar
                                    #install java
                                    apt-get install -y openjdk-8-jre-headless
                                    #run womtool
                                    java -jar womtool-48.jar validate wdl/betamatch_github.wdl -i wdl/betamatch_github.json"""}
        }
      }
    }
    //stage for code metrics etc. Using sonar
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