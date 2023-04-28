import javax.websocket.*;
import javax.websocket.server.ServerEndpoint;
import java.io.IOException;
import java.util.LinkedList;
import java.util.Scanner;
import java.util.ArrayList;

import TN93.*;

@ServerEndpoint("/")
public class Server {
    public Server(){
        System.out.println("class loaded " + this.getClass());
    }

    @OnOpen
    public void onOpen(Session session) {
        System.out.printf("Session opened, id: %s%n", session.getId());
/*        try {
            session.getBasicRemote().sendText("Hi there, we are successfully connected.");
        } catch (IOException ex) {
            ex.printStackTrace();
        }
*/
    }

    @OnMessage
    public void onMessage(String message, Session session) {
        System.out.printf("Message received. Session id: %s Message: %s%n",
                session.getId(), message);
        try {
            ArrayList<Seq> seqs = TN93.read_seqs(new Scanner(message));
            TN93 tn93 = new TN93();
            //TODO: This was intentinally broken when making the disk-memory based TN93
            // Since we no longer return a matrix, rather we write to a file.
            double [][] dist = new double[1][1];

            for (int i = 1; i < dist.length; ++i) {
                for (int j = 0; j < i; ++j) {
                    session.getBasicRemote().sendText(String.format("%s,%s,%f", seqs.get(i).getName(),
                            seqs.get(j).getName(), dist[i][j]));
                }
            }
            //session.getBasicRemote().sendText(String.format("result: %f%n", dist[0][1]));
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

    @OnError
    public void onError(Throwable e) {
        e.printStackTrace();
    }

    @OnClose
    public void onClose(Session session) {
        System.out.printf("Session closed with id: %s%n", session.getId());
    }
}
