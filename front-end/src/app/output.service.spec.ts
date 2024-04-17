import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { DomSanitizer } from '@angular/platform-browser';
import { fakeAsync, tick, TestBed } from '@angular/core/testing';
import { first, firstValueFrom } from 'rxjs';
import { HttpClientTestingModule, HttpTestingController } from '@angular/common/http/testing';
import { MatSnackBar, MatSnackBarModule } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { MockUserIdService } from './mock-user-id.service';
import { OutputService } from './output.service';
import { UserIdService } from './user-id.service';

import { createServer, Server } from "http";
import { AddressInfo } from "net";
import { io as ioc, Socket } from "socket.io-client";
import { BackendSocketClient } from './backend-socket.client';
import { MockBackendSocketClient } from './mock-backend-socket.client';

describe('OutputService', () => {
  let service: OutputService;
  let snackBar: MatSnackBar;
  let sanitizer: DomSanitizer;
  let io: Server, serverSocket: Socket;

  beforeEach(() => {
    TestBed.configureTestingModule({
      imports: [
        BrowserAnimationsModule,
        HttpClientTestingModule,
        MatSnackBarModule,
      ], 
      providers: [
        { provide: UserIdService, useClass: MockUserIdService },
        { provide: BackendSocketClient, useClass: MockBackendSocketClient }
      ],
    });
    service = TestBed.inject(OutputService);
    snackBar = TestBed.inject(MatSnackBar);
    sanitizer = TestBed.inject(DomSanitizer);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });

  describe('outputs', () => {
    it('should initially be empty', async () => {
      const outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('resetOutputs', () => {
    it('should set the outputs list to be empty', async () => {
      let outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(1);
      service.resetOutputs();
      outputs = await firstValueFrom(service.outputs.pipe(first()));
      expect(outputs.length).toBe(0);
    });
  });

  describe('executeBlocks', () => {
    it('should open snack bar when userId is null', () => {
      const spy = spyOn(snackBar, 'open');
      const blocks: Block[] = [];
      service.executeBlocks(blocks);
      expect(spy).toHaveBeenCalledOnceWith('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    });

    it('should call backendSocketClient.sendRequest() with the appropriate message', () => {
      const userIdService: UserIdService = TestBed.inject(UserIdService);
      userIdService.setUserId();
      const backendSocketClient: BackendSocketClient = TestBed.inject(BackendSocketClient);
      const spy = spyOn(backendSocketClient, 'sendRequest');
      const block: Block = {
        blockId: 'loaddata',
        title: 'Load Data',
        possibleChildBlocks: [],
        parameters: [
          {key: 'foo', text: 'bar', value: 42},
        ],
      };
      const blocks: Block[] = [];
      blocks.push(block);
      service.executeBlocks(blocks);
      const message = {
        user_id: 'mock_user_id',
        blocks: [{ block_id: 'loaddata', foo: 42}],
      };
      expect(spy).toHaveBeenCalledOnceWith(message);
    });



    it('should add a response to outputs array when it receives a valid response', fakeAsync(() => {
      // const block: Block = {
      //   blockId: 'loaddata',
      //   title: 'Load Data',
      //   possibleChildBlocks: [],
      //   parameters: [],
      // };
      // const userIdService: UserIdService = TestBed.inject(UserIdService);
      // userIdService.setUserId();
      // const mockHttp = TestBed.inject(HttpTestingController);
      // service.executeBlock(block).then(async () => {
      //   const outputs = await firstValueFrom(service.outputs);
      //   expect(outputs).toEqual([{text: 'Hello World'}]);
      // });
      // tick();
      // const req = mockHttp.expectOne('http://localhost:5000/api/loaddata?user_id=mock_user_id');
      // expect(req.request.method).toBe('GET');
      // req.flush({text: 'Hello World'});
      // mockHttp.verify();
      expect(false).toBe(true);
    }));

    it('should replace image outputs with a sanitised SafeUrl', fakeAsync(() => {
      // const block: Block = {
      //   blockId: 'qcplots',
      //   title: 'Quality Control Plots',
      //   possibleChildBlocks: [],
      //   parameters: [],
      // };
      // const userIdService: UserIdService = TestBed.inject(UserIdService);
      // userIdService.setUserId();
      // const mockHttp = TestBed.inject(HttpTestingController);
      // const spy = spyOn(sanitizer, 'bypassSecurityTrustUrl');
      // service.executeBlock(block).then(async () => {
      //   expect(spy).toHaveBeenCalledTimes(1);
      //   const outputs = await firstValueFrom(service.outputs);
      //   expect(outputs.length).toBe(1);
      //   const expectedValue = sanitizer.bypassSecurityTrustUrl('data:image/png;base64,' + 'Image string');
      //   expect(outputs[0].img).toBe(expectedValue);
      //   expect(outputs[0].alttext).toBe('Alt text');
      // });
      // tick();
      // const req = mockHttp.expectOne('http://localhost:5000/api/qcplots?user_id=mock_user_id');
      // expect(req.request.method).toBe('GET');
      // req.flush({img: 'Image string', alttext: 'Alt text'});
      // mockHttp.verify();
      expect(false).toBe(true);
    }));

    it('should open snack bar with appropriate message when response is a 406 error', fakeAsync(() => {
      // const block: Block = {
      //   blockId: 'loaddata',
      //   title: 'Load Data',
      //   possibleChildBlocks: [],
      //   parameters: [],
      // };
      // const userIdService: UserIdService = TestBed.inject(UserIdService);
      // userIdService.setUserId();
      // const mockHttp = TestBed.inject(HttpTestingController);
      // const spy = spyOn(snackBar, 'open');
      // service.executeBlock(block).then(async () => {
      //   expect(spy).toHaveBeenCalledOnceWith('The blocks you have executed are not a valid order. Please check the blocks and try again.', 'Close', { duration: 5000 });
      // });
      // tick();
      // const req = mockHttp.expectOne('http://localhost:5000/api/loaddata?user_id=mock_user_id');
      // expect(req.request.method).toBe('GET');
      // req.flush('', { status: 406, statusText: 'Bad Request'});
      // mockHttp.verify();
      expect(false).toBe(true);
    }));

    it('should open snack bar with appropriate message for other errors', fakeAsync(() => {
      // const block: Block = {
      //   blockId: 'runumap',
      //   title: 'Run UMAP',
      //   possibleChildBlocks: [],
      //   parameters: [],
      // };
      // const userIdService: UserIdService = TestBed.inject(UserIdService);
      // userIdService.setUserId();
      // const mockHttp = TestBed.inject(HttpTestingController);
      // const spy = spyOn(snackBar, 'open');
      // service.executeBlock(block).then(async () => {
      //   expect(spy).toHaveBeenCalledOnceWith('There has been an error. Please refresh the page and try again.', 'Close', { duration: 5000 });
      // });
      // tick();
      // const req = mockHttp.expectOne('http://localhost:5000/api/runumap?user_id=mock_user_id');
      // expect(req.request.method).toBe('GET');
      // req.flush('', { status: 400, statusText: 'Bad Request'});
      // mockHttp.verify();
      expect(false).toBe(true);
    }));
  }); 
});
